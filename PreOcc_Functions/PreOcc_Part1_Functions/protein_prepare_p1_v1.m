function [median_vals,sum_vals,out_filename,peps_keep_inds] = ...
    protein_prepare_p1_v1(filename_in,measurement_chans, num_channels,date_out)
    
unmod_peptides = readtable(filename_in,'ReadVariableNames',1);
%unmod_peptides = readtable(filename_in,'delimiter','\t','headerlines',0,'ReadVariableNames',1);

find_fileinfo = split(filename_in,'.csv');
fileinfo = find_fileinfo{1};

sum_sn_cutoff = 0; 
ref_name = {'ProteinId'};
pep_seq_name = {'PeptideSequence'};
lin_2_inds = 1:height(unmod_peptides);

unmod_unfiltered_sn = unmod_peptides{:,measurement_chans(1:num_channels)};
sum_sn = sum(unmod_unfiltered_sn,2);
peps_keep_sn = logical(sum_sn >= sum_sn_cutoff);

unmod_refs = unmod_peptides{:,ref_name};

peps_contams = contains(unmod_refs,'contam');
peps_not_contams = ~peps_contams;

peps_fdr = strncmp(unmod_refs,'##',2);
peps_not_fdr = ~peps_fdr;

% pep_seq = unmod_peptides{:,pep_seq_name};
% peps_oxMet = contains(pep_seq,'M*');
% peps_not_oxMet = ~peps_oxMet;

peps_keep = all([peps_keep_sn peps_not_contams peps_not_fdr],2);

peps_keep_inds = lin_2_inds(peps_keep);

unmod_peps_filter = unmod_peptides(peps_keep,:);
%refs_keep = unmod_refs(peps_keep);

% out_pep_GS = match_refs2syms(refs_keep,all_refs);
% pep_GS = cell(length(out_pep_GS),1);
% pep_GS(logical(out_pep_GS)) = all_syms(out_pep_GS(logical(out_pep_GS)));
% unmod_peps_filter.Gene_Symbol = pep_GS;

clear unmod_peptides

% now compute the median norm of the peptides that you are keeping 
median_vals = median(unmod_peps_filter{:,measurement_chans(1:num_channels)},1);
sum_vals = sum(unmod_peps_filter{:,measurement_chans(1:num_channels)},1); 
% 
% now we want to find the collumn that has the protein reference
% information and replace it with Protein_Reference 

ind_ref = find(strcmp(ref_name,unmod_peps_filter.Properties.VariableNames));
unmod_peps_filter.Properties.VariableNames{ind_ref} = 'Protein_Reference';

% now we want to find the column that has the peptide sequence 
% information and replace it with Peptide_Sequence

ind2_ref = find(strcmp(pep_seq_name,unmod_peps_filter.Properties.VariableNames));
unmod_peps_filter.Properties.VariableNames{ind2_ref} = 'Peptide_Sequence';

% now we want to find the column that has the phospho site position on the
% reference and save it as Phosphosite_Position

data_col_num = num_channels;
Channel_Names = cell([1,data_col_num]);

for i = 1:data_col_num
    New_Channel_Name = strcat('Chan_',num2str(i));
    unmod_peps_filter.Properties.VariableNames{measurement_chans(i)} = ...
    New_Channel_Name;
    Channel_Names{i} = New_Channel_Name;
end

clear unmod_peptides

unmod_data = unmod_peps_filter{:,Channel_Names};
unmod_data_r = round(unmod_data,1);

unmod_peps_info = unmod_peps_filter(:,{'Protein_Reference','Peptide_Sequence'});

unmod_peps_data = array2table(unmod_data_r,'VariableNames',Channel_Names);

out_uT = [unmod_peps_info unmod_peps_data];
out_uT_S = sortrows(out_uT,[1 2]);


out_filename = [pwd,filesep,fileinfo,'_forMatching','.csv'];
%writetable([unmod_peps_info unmod_peps_data],out_filename);
writetable(out_uT_S,out_filename);

end