function [out_filename_keep,num_peps_REM] = ...
    PreOcc_Major_One_Phos(input_file,filename_reference,...
    Channel_Names,median_scale,dir_str,exp_str)

%   mod_peps_info_REMOVE,mod_data_REMOVE,phos_median

%  INPUT FILE has the following columns
% Protein_Reference, Peptide_Sequence, Phosphosite_Position,
% Phosphosite_Motif, Phosphosite_LocalizeScore, GeneSymbol 
% and then potentially other columsn with GeneSymbol information

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
ModPepImport = readtable(input_file,'format','auto');
initial_inds_mod = 1:height(ModPepImport);

% SPECIFICALLY load the protein references, the peptide sequence, and the
% phosphosite position from the input file 

phos_protein = ModPepImport{:,'Protein_Reference'};
phos_position = ModPepImport{:,'Phosphosite_Position'};
phos_sequence = ModPepImport{:,'Peptide_Sequence'};

%%% NOW WE ARE GOING TO GO THROUGH all the REFERENCES and then all the
%%% phosphosite on that reference 

[unique_prot_refs,ia_PR,ic_PR] = unique(phos_protein,'stable');
store_keep_PHOS = zeros(height(ModPepImport),1);
store_phos_POSREL = cell(height(ModPepImport),1);

prov_PhosForms_counter = 1;
set_counter = 1;

set_indexing_array = [];

for i = 1:length(unique_prot_refs)
    
    phos_protIdindex = strcmp(phos_protein(ia_PR(i)),FASTA_proteinId);
    ProteinSequence = FASTA_sequence(phos_protIdindex);
    ProteinName = phos_protein(ia_PR(i));
    
    phos_position_PR = phos_position(ic_PR==i); % find the positions of all of 
    % phospho-peptides assigned to this reference 
    phos_PEPS_prot = phos_sequence(ic_PR == i);
    
    [unique_positionS,ia_PP,ic_PP] = unique(phos_position_PR,'stable');
    inds_2_store = initial_inds_mod(ic_PR==i);
    store_span = zeros(5,length(unique_positionS));
    
    for j = 1:length(unique_positionS) 
    % first remove the miss-cleaved peptides 
    
        phos_peps_1 = phos_PEPS_prot(ic_PP==j);
        store_span(1,j) = prov_PhosForms_counter;
        
        if (size(phos_peps_1,1)>1) % see if there are differential missed cleavages
            
           % phos_peps_2 = regexprep(phos_peps_1,'#|*|^[-|R|K].|.[A-Z|-]$','');
              phos_peps_2 = regexprep(phos_peps_1,'@|#|*|^[-|R|K].|.[A-Z|-]$','');
            [uu_ph,~,uu_ic] = unique(phos_peps_2);
            store_rk_ph = zeros(length(phos_peps_2),1);
            
            for k = 1:length(uu_ph)
                check_p = uu_ph{k,:};
                check_rk = regexp(check_p(1:end-1),'[R|K]','start');
                if (isempty(check_rk))
                else
                    store_rk_ph(uu_ic==k) = length(check_rk);
                end
            end
            
            if (min(store_rk_ph)==max(store_rk_ph)) % if all peptides have same num MC
                store_keep_PHOS(inds_2_store(ic_PP==j)) = prov_PhosForms_counter;
            else
                min_v = min(store_rk_ph);
                ind_all_IC = find(ic_PP==j);
                ind_all_min = ind_all_IC(store_rk_ph==min_v);
                store_keep_PHOS(inds_2_store(ind_all_min)) = prov_PhosForms_counter;
            end
            
        else
            store_keep_PHOS(inds_2_store(ic_PP==j)) = prov_PhosForms_counter; % keep the single phospho-measuremen
        end
        
        % use the sequence of one of the kept peptides to find the longest
        % spanning sequence 
        ind_SEQ = find(store_keep_PHOS == prov_PhosForms_counter,1);
        %phos_seq = regexprep(phos_sequence(ind_SEQ),'#|*|^[R|K].|.[A-Z|-]$','');
        phos_seq = regexprep(phos_sequence(ind_SEQ),'@|#|*|^[R|K].|.[A-Z|-]$','');
        peplength = length(char(phos_seq));

        %Search peptide against protein sequence to find start position
        pepstart = regexp(ProteinSequence,phos_seq);
        
        %I and L are sometimes substituted...if something does match this is
        %likely why...check for this
        if isempty(pepstart)
            
            if isempty(regexp(phos_seq,'I|L','ONCE'))==0
                pepchar_IL = regexprep(phos_seq,'I|L','[IL]');
                pepstart = regexp(ProteinSequence,pepchar_IL);
                
                if (isempty(pepstart))
                     warning(['Protein Reference Flag- Peptide', char(phos_seq) '" failed to align to the sequence from ProteinId "',char(ProteinName),'". This usually indicates a mismatch between the protein reference provided as input and the reference that was used to search the mass spectrometry data. Check to see if the correct input was used. It can also mean there is a formating mismatch between the ProteinId in the fasta file versus the ouput MS data.'])
                     pepstart = 0;
                     peplength = 0;
                else
                    pepstart = double(pepstart{1,1});
                end
            else
                %warn1 = UnmodSpecificProtein{peptidecounter,NonData_VarNames{1}}{1,1};
                warning(['Protein Reference Flag- Peptide "', phos_seq '" failed to align to the sequence from ProteinId "',warn1,'". This usually indicates a mismatch between the protein reference provided as input and the reference that was used to search the mass spectrometry data. Check to see if the correct input was used. It can also mean there is a formating mismatch between the ProteinId in the fasta file versus the ouput MS data.'])
                pepstart = 0;
                peplength = 0;
            end
            
        else
            pepstart = double(pepstart{1,1});
            if (size(pepstart,2)>1)
                disp(['oh no, this peptide aligns to many locations on protein ',char(ProteinName),'.'])
                pepstart = pepstart(1);
            end
        end
        
        if (isempty(pepstart))
            
        else
            store_span(2,j) = pepstart;
            store_span(3,j) = store_span(2,j)+peplength-1;
            
           % find the position of the phosphrylated residue relative to the
           % peptide 
           
           if (store_span(2,j)>0)
               pos_info = unique_positionS{j};
               pos_info_array = split(pos_info,';');
               pos_update = '';
               
               for M = 1:length(pos_info_array)
                   
                   pos_num = str2double(pos_info_array{M});
                   pos_num_delta = pos_num - store_span(2,j)+1; % find how many it is after the start position
                   
                   
                   if (pos_num_delta > store_span(3,j))
                       disp('Oh! Looks like removing MC peptides resulted in loss of this phospho-residue')
                   else
                       pos_update = [pos_update,';',num2str(pos_num_delta)];
                       
                   end
                   
               end
           else
               pos_update = '';
               
           end
            
            
        end

%             %find range of residue positions 
%             UnmodResiduePositions(peptidecounter,1:length(pepchar{1,1})) = ...
%                 pepstartNum:1:((pepstartNum-1)+peplength);    
        store_phos_POSREL{prov_PhosForms_counter} = pos_update;
        prov_PhosForms_counter = prov_PhosForms_counter + 1;
        
        
    end
    
    store_span_S = sortrows(store_span',2);
    % sort store_span_S based on the start position of the sequence 
        
    for k = 1:size(store_span_S,1)
        set_up_flag = 0;
        
        s_check = store_span_S(k,2);
        s_lin_inds = k+1:1:size(store_span_S,1);
        s_other_inds = logical(store_span_S(s_lin_inds,4));
        s_others = store_span_S(s_lin_inds(~s_other_inds),:);
        
        if logical(store_span_S(k,4))
            % then we have already assigned a set number
            
        elseif (isempty(s_others))
            store_span_S(k,4) = set_counter;
            set_counter = set_counter + 1;
            
        else
            
            out_EQ_S = find(s_check == s_others(:,2));
            
            e_check = store_span_S(k,3);
            out_EQ_E = find(e_check == s_others(:,3));
            out_OL = find(e_check > s_others(:,2));
            
            if (out_EQ_S) % if two or more peptides have same starting residue
                
                store_span_S(k,4) = set_counter;
                store_span_S(out_EQ_S+k*ones(size(out_EQ_S,1),1),4) = set_counter;
                set_up_flag = 1;
%                set_counter = set_counter + 1;
            end 
            
            if (out_EQ_E) % if two or more peptides have the same ending residue 
                
                store_span_S(k,4) = set_counter;
                store_span_S(out_EQ_E+k*ones(size(out_EQ_E,1),1),4) = set_counter;
                set_up_flag = 1;
            end    
                
            if (out_OL) % if the last positoin of the peptide of interest is larger than the
                % starting residue of other peptides then those peptides
                % overlap
                
                store_span_S(k,4) = set_counter;
                store_span_S(out_OL+k*ones(size(out_OL,1),1),4) = set_counter;
                set_up_flag = 1;
            end
            
            if (set_up_flag == 1)
                set_counter = set_counter +1;
            else
                store_span_S(k,4) = set_counter;
                set_counter = set_counter + 1;
            end
            
        end
            
       
    end
    
    store_span_S_sort = sortrows(store_span_S,4);
    set_indexing_array = [set_indexing_array; store_span_S_sort(:,[1 4])]; 
    
end

% first find the information that we are going to keep and keep and sort it

phos_keep_inds = (logical(store_keep_PHOS));

store_p_T = array2table(store_keep_PHOS,'VariableNames',{'PhosForm_1'});
ModPepImport_Save = [store_p_T ModPepImport];
ModPepImport_Save = ModPepImport_Save(phos_keep_inds,:);
ModPepImport_Save = sortrows(ModPepImport_Save,1);

ModPepImport_REM = ModPepImport(~phos_keep_inds,:);
num_peps_REM = sum(~phos_keep_inds);

% now we need to through and assign a set #, phospho-form IDENTIFIER, and
% info for where PHOSPHORYLATION is relative to the peptide sequence 
out_SETS = zeros(height(ModPepImport_Save),1);
out_FORMS = cell(height(ModPepImport_Save),3);

[u_forms,ia_forms,ic_forms] = unique(ModPepImport_Save{:,1});
peps_out_counter = 1;
num_max_ints = length(char(num2str(length(u_forms))));
spacer = cell(1,num_max_ints);
spacer(:) = {'_'};

for i = 1:length(u_forms)
    
    p_form = u_forms(i); % this is the phospho-form of interest
    si_ind = find(set_indexing_array(:,1)==p_form);
    out_SETS(peps_out_counter:peps_out_counter+length(find(ic_forms==i))-1) = ...
        set_indexing_array(si_ind,2);
    
    out_end = cell(1,length(find(ic_forms==i)));
    for j = 1:length(out_end)
        out_end{j} = num2str(p_form);
    end
    num_space = num_max_ints-length(char(num2str(p_form)));
    if (num_space)
        spacer_use = strcat(spacer{1:num_space});
    else
        spacer_use = '';
    end
    out_FORMS(peps_out_counter:peps_out_counter+length(find(ic_forms==i))-1,1) = ...
        strcat('Form_',spacer_use,out_end);
    
    out_FORMS(peps_out_counter:peps_out_counter+length(find(ic_forms==i))-1,2) = ...
        store_phos_POSREL(p_form);  
    
    phos_pep_out = ModPepImport_Save{ia_forms(i),'Peptide_Sequence'};
    phos_pep_out = regexprep(phos_pep_out,'#','');
    phos_pep_out = regexprep(phos_pep_out,'#|@|*','');
    
    out_FORMS(peps_out_counter:peps_out_counter+length(find(ic_forms==i))-1,3) = ...
        phos_pep_out;  
    
    peps_out_counter = peps_out_counter + length(find(ic_forms==i));
end

% scale the data

mod_data_orig = ModPepImport_Save{:,Channel_Names};
mod_data_scaled = mod_data_orig./repmat(median_scale,size(mod_data_orig,1),1); % scale the raw signal to noise values 

mod_data_r = round(mod_data_scaled,1); % round the scaled signal to noise value to 1 
ModPepImport_Save{:,Channel_Names} = mod_data_r;

% scale the removed data  
mod_data_rem = ModPepImport_REM{:,Channel_Names};
mod_data_scaled_rem = mod_data_rem./repmat(median_scale,size(mod_data_rem,1),1); % scale the raw signal to noise values 

mod_data_rem_R = round(mod_data_scaled_rem,1); % round the scaled signal to noise value to 1 
ModPepImport_REM{:,Channel_Names} = mod_data_rem_R;

% round the localize score for easier usability 
localize_phos_score = ModPepImport_Save{:,'Phosphosite_LocalizeScore'};
ModPepImport_Save{:,'Phosphosite_LocalizeScore'} = round(localize_phos_score,1);

% write output file with ModPepImport_Save
% write output file with ModPepImport_REM 
new_vars = [array2table(out_SETS,'VariableNames',{'P_Set'}) ...
    cell2table(out_FORMS,'VariableNames',{'P_Form','PhosPosition_Rel','Peptide_Sequence'})];

old_vars = ModPepImport_Save(:,[{'GeneSymbol','GeneSymbol_Xenbase','Protein_Reference','Phosphosite_Position','Phosphosite_LocalizeScore','Phosphosite_Motif','GeneSymbol_More'},Channel_Names]);

out_T = [new_vars(:,1:2) old_vars(:,1:6) new_vars(:,3:4) old_vars(:,7:end)];
out_T = sortrows(out_T,1);

date_out = datestr(now,'yymmdd'); 

out_filename_keep = [pwd,filesep,dir_str,filesep,date_out,'_',exp_str,'_PHOS_ONE.csv']; % create the name for the new file
writetable(out_T,out_filename_keep); % write combined information and data to the file 

out_filename_remove = [pwd,filesep,dir_str,filesep,date_out,'_',exp_str,'_PHOS_removed.csv']; % create the name for the new file
writetable(ModPepImport_REM,out_filename_remove); % write combined information and data to the file 



end

