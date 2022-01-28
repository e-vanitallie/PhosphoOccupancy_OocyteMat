function [out_filename,Channel_Names] = ...
    phos_prepare_p1_v4(phos_peps_analyze_file,phos_CompositeInfo_search_file,...
    phos_measurement_chans,dir_str,exp_str)
  
jgi_v9p2_genesym = readtable('210113_xenla_v9p2_humangn.xlsx');
jgi_v9p2_refs = jgi_v9p2_genesym{:,'ProteinId'};

Xenbase_names = jgi_v9p2_genesym{:,'laevisGene'};
Human_syms = jgi_v9p2_genesym{:,'humanGene'};
Human_more_info = jgi_v9p2_genesym{:,'description'};

mod_peptides = readtable(phos_peps_analyze_file,'ReadVariableNames',1,'format','auto');


composite_file = readtable(phos_CompositeInfo_search_file,'headerlines',0,'readvariablenames',1,'delimiter',',');

% section 2 -- process the phospho analyze file to go into modsets 

phos_ref_name = {'protein_id'}; %x


phos_pep_sequence = {'sequence'}; %x 
phos_position = {'site_position'}; %x
phos_motif = {'motif_peptide'}; %x
phos_score = {'max_score'}; %x
phos_siteid = {'site_id'}; %x 
compphos_sitepos = {'siteIDstr','sitePosStr','motifPeptideStr'};

phos_refs = mod_peptides{:,phos_ref_name};

% Identify peptides where sun signal to noise is > 200 AND NOT (contam or  ##) 

sum_sn_cutoff = 200; 

phos_unfiltered_sn = mod_peptides{:,phos_measurement_chans};
phos_sum_sn = sum(phos_unfiltered_sn,2);
phos_peps_keep_sn = logical(phos_sum_sn >= sum_sn_cutoff);

phos_peps_contams = contains(phos_refs,'contam');
phos_peps_not_contams = ~phos_peps_contams;

phos_peps_fdr = strncmp(phos_refs,'##',2);
phos_peps_not_fdr = ~phos_peps_fdr;

phos_peps_keep = all([phos_peps_keep_sn phos_peps_not_contams phos_peps_not_fdr],2);

% Cross reference the ProteinIds with Human Gene Symbols 

mod_peps_filter = mod_peptides(phos_peps_keep,:);
mod_refs_keep = phos_refs(phos_peps_keep);

refs_parts = split(mod_refs_keep, '|');
mod_refs_keep_Match = strcat(refs_parts(:,1),'|',refs_parts(:,2),'|',refs_parts(:,4));

% calling another function -> match_refs2syms
out_pep_GS = match_refs2syms(mod_refs_keep_Match,jgi_v9p2_refs);

pep_GS = cell(length(out_pep_GS),1);
pep_GS(logical(out_pep_GS)) = Human_syms(out_pep_GS(logical(out_pep_GS)));

inds_o = zeros(length(pep_GS),1);

for i = 1:length(pep_GS)
    o = find(strcmp(pep_GS{i},Human_syms));
    if(o)
        inds_o(i) = o(1);
    end
end

pep_GS_more = cell(length(inds_o),1);
pep_GS_more(:,1) = {'no_description'}; %initialize the cell so there is always something 
pep_GS_more(logical(inds_o)) = Human_more_info(inds_o(logical(inds_o)));

pep_GS_xenbase = cell(length(out_pep_GS),1);
pep_GS_xenbase(logical(out_pep_GS)) = Xenbase_names(out_pep_GS(logical(out_pep_GS)));

mod_peps_filter.GeneSymbol = pep_GS;
mod_peps_filter.GeneSymbol_More = pep_GS_more;
mod_peps_filter.GeneSymbol_Xenbase = pep_GS_xenbase;

% Update the columns of the output files with standardized names 

ind_ref1 = find(strcmp(phos_siteid,mod_peps_filter.Properties.VariableNames));
mod_peps_filter.Properties.VariableNames{ind_ref1} = 'Site_Id';

% now we want to find the collumn that has the protein reference
% information and replace it with Protein_Reference 

ind_ref = find(strcmp(phos_ref_name,mod_peps_filter.Properties.VariableNames));
mod_peps_filter.Properties.VariableNames{ind_ref} = 'Protein_Reference';

% now we want to find the column that has the peptide sequence 
% information and replace it with Peptide_Sequence

ind2_ref = find(strcmp(phos_pep_sequence,mod_peps_filter.Properties.VariableNames));
mod_peps_filter.Properties.VariableNames{ind2_ref} = 'Peptide_Sequence';

% now we want to find the column that has the phospho site position on the
% reference and save it as Phosphosite_Position

ind3_ref = find(strcmp(phos_position,mod_peps_filter.Properties.VariableNames));
mod_peps_filter.Properties.VariableNames{ind3_ref} = 'Phosphosite_Position';

% now we want to find the column that has the phosphosite motif and save it
% as Phosphosite_Motif

ind4_ref = find(strcmp(phos_motif,mod_peps_filter.Properties.VariableNames));
mod_peps_filter.Properties.VariableNames{ind4_ref} = 'Phosphosite_Motif';

% now we want to find the column that has the phosphosite position score and save it
% as Phosphosite_Localize_Score
ind5_ref = find(strcmp(phos_score,mod_peps_filter.Properties.VariableNames));
mod_peps_filter.Properties.VariableNames{ind5_ref} = 'Phosphosite_LocalizeScore';


%%%% need to combine the site position information in the composite sites
%%%% w/ the information in the main peptide file -- THIS USED TO HAPPEN IN
%%%% modsets_return_inds BUT that is late for something core specific
%%%% within our pipeline 

single_sites_index = cellfun(@isempty,regexp(mod_peps_filter.Site_Id,';'));
    % find the site_ids that DO NOT have the ; and are therefor single
    % sites 
single_site_pos_array = mod_peps_filter{single_sites_index,'Phosphosite_Position'};
    % now we have the positions of ONLY the single sites 
    
single_site_pos_string = cell(size(single_site_pos_array,1),1); 
    % creat a cell to store the single site positions 
for num2str_counter = 1:size(single_site_pos_array,1) 
    single_site_pos_string(num2str_counter,1) = ...
        cellstr(num2str(single_site_pos_array(num2str_counter)));
end

[single_sites_unique, Ia, Ic] = unique(mod_peps_filter{single_sites_index,'Site_Id'},'stable');
%get the unique single site_ids 

single_site_pos_string_unique = single_site_pos_string(Ia); 
% now we have the position for each of the unique single sites 

site2positionComposite = composite_file{:,compphos_sitepos};

%combine this with the equivalent composite site info from the direct
%import
site2position = [[single_sites_unique single_site_pos_string_unique]; [site2positionComposite(:,1:2)]];
multicheck_motif = [zeros(1,length(single_sites_unique)) 1:length(site2positionComposite)];
s2p_Ids = site2position(:,1);
s2p_Position = site2position(:,2);
pp_SiteIds = mod_peps_filter.Site_Id;

%Use the combined site 2 position map to label all the siteids with positions
All_positions = cell(size(pp_SiteIds,1),1);
All_motifs = mod_peps_filter{:,'Phosphosite_Motif'};

for index = 1:size(pp_SiteIds,1)
    
    a = find(strcmp(pp_SiteIds(index),s2p_Ids));
    
    if (a)
        All_positions(index,1) = s2p_Position(a);
        if (multicheck_motif(a))
            All_motifs(index) = site2positionComposite(multicheck_motif(a),3);
        end
    elseif isempty(a)
        disp('oh no!, we cant find position information for this Site_Id')
        
    end
end

mod_peps_filter.Phosphosite_Position = All_positions; 
mod_peps_filter.Phosphosite_Motif = All_motifs;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now finished with that section


data_col_num = length(phos_measurement_chans);
Channel_Names = cell([1,data_col_num]);

for i = 1:data_col_num
    New_Channel_Name = strcat('Chan_',num2str(i));
    mod_peps_filter.Properties.VariableNames{phos_measurement_chans(i)} = ...
    New_Channel_Name;
    Channel_Names{i} = New_Channel_Name;
    
end

mod_data = mod_peps_filter{:,Channel_Names};

mod_peps_info = mod_peps_filter(:,{'Protein_Reference','Peptide_Sequence','Phosphosite_Position',....
    'GeneSymbol','GeneSymbol_Xenbase','GeneSymbol_More','Phosphosite_Motif','Phosphosite_LocalizeScore'});

[mod_peps_info,s_indexes] = sortrows(mod_peps_info,1);
mod_data = mod_data(s_indexes,:);

mod_peps_data_T = array2table(mod_data,'VariableNames',Channel_Names);

date_out = datestr(now,'yymmdd'); 
out_filename = [pwd,'/',dir_str,'/',date_out,'_',exp_str,'_PHOS_NO_CORE.csv']; % create the name for the new file
writetable([mod_peps_info mod_peps_data_T],out_filename); % write combined information and data to the file 

end

