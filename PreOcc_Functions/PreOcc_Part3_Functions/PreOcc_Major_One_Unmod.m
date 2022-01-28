function [protein_trend_T,parent_trend_T,REMOVED_T,store_matched_sets,mc_peps_rem] =   ...
    PreOcc_Major_One_Unmod(filename_Mod,filename_Unmod,...
    filename_reference,median_scale_in,data_names,chans_info)

% LOAD the fasta reference
PHOS_fasta = fastaread(filename_reference);

FASTA_import = struct2table(PHOS_fasta);
FASTA_import.Properties.VariableNames = {'Header','Sequence'};

FASTA_sequence = FASTA_import.Sequence;
FASTA_Header = FASTA_import.Header;

%Some Id entries have additional information other than just the id, delimit
% on the space to get rid of the superfluous characters
splitid = regexp(FASTA_Header,' ','split');

FASTA_proteinId = cell(size(splitid,1),1);
for Retrieve_splitId = 1:size(splitid)
    FASTA_proteinId(Retrieve_splitId,1) = splitid{Retrieve_splitId,1}(1);
end

% INITIALIZE for STORING
UnmodPepImport = readtable(filename_Unmod,'format','auto');

% SPECIFICALLY load the protein references, the peptide sequence from the input file

pa_protein = UnmodPepImport{:,'Protein_Reference'};
pa_sequence = UnmodPepImport{:,'Peptide_Sequence'};

% determine the vector of unique protein references and have the indexes
%   to peptides associated with those references

[unique_prot_refs_unmod,~,ic_PA] = unique(pa_protein,'stable');
store_parent_UNMOD = zeros(height(UnmodPepImport),1);
% will have the parent set number in the position if it is assigned to
% parent set
store_parent_PROT = ones(height(UnmodPepImport),1);
% will change to zero if the peptide should not be used for protein
% de-trending

%%% NOW WE ARE GOING TO GO THROUGH all the REFERENCES and then all the
%%% phospho-sets on that reference to look for parent sets
SetImport = readtable(filename_Mod,'format','auto');
ph_protein = SetImport{:,'Protein_Reference'};
ph_position = SetImport{:,'Phosphosite_Position'};
ph_sets = SetImport{:,'P_Set'};
[unique_prot_refs_mod,~,ic_PH] = unique(ph_protein,'stable');
store_matched_sets = zeros(max(ph_sets),1);
% 0 if set is unmathced, 1 if set is matched

for i = 1:length(unique_prot_refs_mod)
    
    ph_ref = unique_prot_refs_mod(i);
    sets_ref = ph_sets(ic_PH==i);
    [sets_u,~,ic_sets] = unique(sets_ref);
    ph_positions_ref = ph_position(ic_PH==i);
    
    unmod_ind = find(strcmp(ph_ref,unique_prot_refs_unmod));
    
    if (isempty(unmod_ind))
    else % there are unmodified peptides associated with this phosphorylation
        
        unmod_protIdindex = strcmp(ph_ref,FASTA_proteinId);
        ProteinSequence = FASTA_sequence(unmod_protIdindex);
        
        inds_U = find(ic_PA==unmod_ind);
        unmod_peps = pa_sequence(inds_U);
%        unmod_peps = regexprep(unmod_peps,'^[-|R|K].|.[A-Z|-]$','');
         unmod_peps = regexprep(unmod_peps,'*|#|^[-|R|K].|.[A-Z|-]$','');
        % removing the flanking amino acids
        store_prov_sets = zeros(4,size(unmod_peps,1));
        
        for j = 1:length(unmod_peps)
            
            store_prov_sets(1,j) = j;
            unmod_seq = unmod_peps{j};
            peplength = length(char(unmod_seq));
            
            %Search peptide against protein sequence to find start position
            pepstart = regexp(ProteinSequence,unmod_seq);
            
            %I and L are sometimes substituted...if something does match this is
            %likely why...check for this
            if isempty(pepstart)
                
                if isempty(regexp(unmod_seq,'I|L','ONCE'))==0
                    pepchar_IL = regexprep(unmod_seq,'I|L','[IL]');
                    pepstart = regexp(ProteinSequence,pepchar_IL);
                    
                    if (isempty(pepstart))
                        warning(['Protein Reference Flag- Peptide', char(unmod_seq) '" failed to align to the sequence from ProteinId "',char(ph_ref),'". This usually indicates a mismatch between the protein reference provided as input and the reference that was used to search the mass spectrometry data. Check to see if the correct input was used. It can also mean there is a formating mismatch between the ProteinId in the fasta file versus the ouput MS data.'])
                        pepstart = 0;
                        peplength = 0;
                    else
                        pepstart = double(pepstart{1,1});
                    end
                else
                    warning(['Protein Reference Flag- Peptide "', char(unmod_seq) '" failed to align to the sequence from ProteinId "',char(ph_ref),'". This usually indicates a mismatch between the protein reference provided as input and the reference that was used to search the mass spectrometry data. Check to see if the correct input was used. It can also mean there is a formating mismatch between the ProteinId in the fasta file versus the ouput MS data.'])
                    pepstart = 0;
                    peplength = 0;
                end
                
            else
                pepstart = double(pepstart{1,1});
                if (size(pepstart,2)>1)
                    disp(['oh no, this peptide ',char(unmod_seq),' aligns to many locations on protein ',char(ph_ref),'.'])
                    pepstart = 0;
                    peplength = 0;
                end
   
            end
            
            if (isempty(pepstart))
                disp(['Oh no, issue with peptide ',char(unmod_seq),' alignment to protein'])
            else
             store_prov_sets(2,j) = pepstart;
             store_prov_sets(3,j) = pepstart+peplength-1;
            end
            
        end
        % remove the peptides that map to multiple places
        inds_rem = logical(store_prov_sets(2,:)==0);
        store_parent_PROT(inds_U(inds_rem)) = 0;
        inds_U = inds_U(~inds_rem);
        
        % sort the peptides by starting position
        [store_prov_sets_S,prov_ind] = sortrows(store_prov_sets(:,~inds_rem)',2);
        inds_U = inds_U(prov_ind);
        store_prov_sets_S = store_prov_sets_S';
        
        % find overlap w/ phospho-residues on sets
        flag = 0;
        k = 1;
        while( k<=length(sets_u) && flag == 0)
            ph_positions_set = ph_positions_ref(ic_sets == k);
            
            % now make a numeric array of the unique phosphorylated
            % residues
            out_nums = [];
            for m = 1:length(ph_positions_set)
                o = split(ph_positions_set(m),';');
                for n = 1:length(o)
                    if (ismember(str2double(o(n)),out_nums))
                    else
                        out_nums = [out_nums str2double(o(n))];
                    end
                end
            end
            inds_check = zeros(1,size(store_prov_sets_S,2));
            
            for p = 1:length(out_nums)
                inds_S = logical(out_nums(p)>= store_prov_sets_S(2,:));
                inds_E = logical(out_nums(p)<= store_prov_sets_S(3,:));
                inds_c = and(inds_S,inds_E);
                inds_check = inds_check + inds_c;
            end
            
            % any peptides of inds_check that are greater than 0 will be
            % removed
            % any peptides that have both the same value as length of the
            % numbers will be considered provional set
            
            p_sets = logical(inds_check == length(out_nums));
            if (sum(p_sets)) % if there are 1s in p_sets
                store_parent_UNMOD(inds_U(p_sets))=sets_u(k);
                store_matched_sets(sets_u(k)) = 1;
            end
            
            % now remove from protein all peptides with that overlap all
            % or some of the residues 
            p_rem = logical(inds_check);
            p_rem = or(p_sets,p_rem);
            store_parent_PROT(inds_U(p_rem)) = 0;
            inds_U = inds_U(~p_rem);
            store_prov_sets_S = store_prov_sets_S(:,~p_rem);
            k = k+1;
            
            if (isempty(inds_U))
                flag = 1;
            end
                
        end
    end
end

% go through all the peptides that have been provsionally assigned to sets
% and remove the miss-cleaved ones 
sets_matched = find(store_matched_sets);
mc_peps_rem = 0;

for i = 1:length(sets_matched)
    
    set_num = sets_matched(i);
    prov_inds_seq = find(store_parent_UNMOD==set_num);
    prov_pa_peps_1 = pa_sequence(prov_inds_seq);
    
    if (size(prov_pa_peps_1,1)>1) % see if there are differential missed cleavages
        
       %prov_pa_peps_2 = regexprep(prov_pa_peps_1,'*|^[-|R|K].|.[A-Z|-]$','');
        prov_pa_peps_2 = regexprep(prov_pa_peps_1,'#|*|^[-|R|K].|.[A-Z|-]$','');
        [uu_pa,~,uu_ic] = unique(prov_pa_peps_2);
        store_rk_pa = zeros(length(prov_pa_peps_2),1);
        
        for k = 1:length(uu_pa)
            check_p = uu_pa{k,:};
            check_rk = regexp(check_p(1:end-1),'[R|K]','start');
            if (isempty(check_rk))
            else
                store_rk_pa(uu_ic==k) = length(check_rk);
            end
        end
        
        if (min(store_rk_pa)==max(store_rk_pa)) % if all peptides have same num MC
            % don't do anything because we are keeping all provisional
            % peptides
        else
            min_v = min(store_rk_pa);
            ind_all_max = prov_inds_seq(~store_rk_pa==min_v);
            store_parent_UNMOD(ind_all_max) = 0;
            mc_peps_rem = mc_peps_rem + length(ind_all_max);
        end
    end
end

% SCALE the data
unmod_data_orig = UnmodPepImport{:,data_names};
unmod_data_scaled = unmod_data_orig./repmat(median_scale_in,size(unmod_data_orig,1),1); % scale the raw signal to noise values 
unmod_data_r = round(unmod_data_scaled,1); % round the scaled signal to noise value to 1 

% create output table with the SCALED pepetides for PROTEIN trend estimation
protein_trend_T = [UnmodPepImport(logical(store_parent_PROT),1:chans_info) array2table(unmod_data_r(logical(store_parent_PROT),:),'VariableNames',data_names)];

% create an output table the SCALED peptides for PARENT trend estimattion 
baciq_pre_unmod = 'Set_';
sets_num_peps = store_parent_UNMOD(logical(store_parent_UNMOD));
sets_array = cell(length(sets_num_peps),1);

num_max_ints = length(char(num2str(max(sets_num_peps))));
spacer = cell(1,num_max_ints);
spacer(:) = {'_'};

for i = 1:length(sets_array)    
    num_space = num_max_ints-length(char(num2str(sets_num_peps(i))));
    if (num_space)
        spacer_use = strcat(spacer{1:num_space});
    else
        spacer_use = '';
    end
    sets_array(i) = cellstr(strcat(baciq_pre_unmod,spacer_use,num2str(sets_num_peps(i))));
end

parent_trend_T = [cell2table(sets_array,'VariableName',{'Set'}) UnmodPepImport(logical(store_parent_UNMOD),1:chans_info) array2table(unmod_data_r(logical(store_parent_UNMOD),:),'VariableNames',data_names)];

% write an output table with the SCALED peptides for assignment to NEITHER 
completely_removed = ~(or(logical(store_parent_PROT),logical(store_parent_UNMOD)));
REMOVED_T = [UnmodPepImport(completely_removed,1:chans_info) array2table(unmod_data_r(completely_removed,:),'VariableNames',data_names)];


end
