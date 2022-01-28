function [file_4_occ,data_store] = ...
    PreOcc_Major_Two_part2...
    (set_info,phos_info,occ_phos_data,occ_parent_data,...
    prot_proteinIDs,prot_data,occ_table_var_types_extra)

% THIS FUNCTION combines the info & data about phospho-sets, matched
% parents, and protein trends into the correct fromat for input into modocc
%
%%%% INPUTS %%%
%
% set_info - # of sets x 5 cell array
%     has te Set #s where have parent peptide data and then the
%     Protein_Reference, GeneSymobl, and GeneSymobl_Xenbase information for
%     the set -- (this information is not different for different
%     phosphoforms which is why it is the set info)
%
% phos_info - # of unique phospho-site_ids that are part of sets w/
% occupancy x 3 cell array
%     has the Set #s to connect the data to Sets, the PhosphoSite_Motif for
%     the specific phospho_site_id and the phophosite_position on the
%     protein reference
%
% occ_phos_data - MATRIX that is size(# of unique phospho-site_ids, # of
% data_columns has the data for creating the input into mod_occ
%
% occ_parent_data - MATRIX that is size(# of matched sets, # of
% data_columns) that has the PARENT data for creating the input into mod_occ
%
% prot_proteinIDs - the cell array of the protein sequence reference
%     identifiers in the same order as the BACIQd data so that we can find
%     the protein data for normalization
%
% prot_data - the matrix of BACIQd non-Parent protein data

% creat a matrix of the correct size and two tables of the correct size

h_newtable = 2*size(occ_parent_data,1)+size(occ_phos_data,1);
file_4_occ = table('Size',[h_newtable 1+width(set_info)],'VariableTypes', ...
    [{'double','string','double','string','string'},occ_table_var_types_extra]);

data_store = []; % zeros(h_newtable,width(occ_phos_data));

occ_phos_inds = phos_info{:,1};
% this tells us the sets that phos occ data is associated with

counter = 1;
data_col_num = size(occ_parent_data,2);

for i = 1:size(occ_parent_data,1)
    
    set_num = set_info{i,1};
    s_prot = set_info{i,'Protein_Reference'};
    
    find_prot_measure = find(contains(prot_proteinIDs,s_prot));
    
    if (find_prot_measure)
        match_prot_log = 1;
        matched_prot_data = prot_data(find_prot_measure(1),:);
      %  prot_u_peps = numUnique_peps_prot(find_prot_measure(1));
        
    else
        match_prot_log = 0;
       % prot_u_peps = 0;
    end
    
    file_4_occ{counter,1} = set_num;
    file_4_occ(counter,2) = {strcat(num2str(set_num),'_unmodparent')};
    file_4_occ{counter,3} = match_prot_log;
    file_4_occ(counter,4) = {s_prot};
    file_4_occ(counter,5) = set_info(i,'GeneSymbol');
    file_4_occ(counter,6) = set_info(i,5);
    file_4_occ(counter,9) = set_info(i,'PP_Reference');
        
    data_store = [data_store; occ_parent_data(i,:)];
    
    % find indexes to associatd phos data
    inds_2_phos = find(occ_phos_inds == set_num);
    
    for j = 1:length(inds_2_phos)
        
        file_4_occ{counter+j,1} = set_num;
        file_4_occ(counter+j,2) = phos_info(inds_2_phos(j),'P_Form');
        file_4_occ{counter+j,3} = match_prot_log;
        file_4_occ(counter+j,4) = {s_prot};
        file_4_occ(counter+j,5) = set_info{i,'GeneSymbol'};
         file_4_occ(counter+j,6) = set_info(i,5);
          file_4_occ(counter+j,7:8) = phos_info(inds_2_phos(j),{'Phosphosite_Position','Phosphosite_Motif'});
          file_4_occ(counter+j,9:11) = phos_info(inds_2_phos(j),{'PP_Reference','Human_Residue','Human_Motif'});
          file_4_occ(counter+j,12) = phos_info(inds_2_phos(j),{'MoreInfoFlag'});
        
        data_store = [data_store; occ_phos_data(inds_2_phos(j),:)];
        
    end
    
    counter = counter+length(inds_2_phos)+1;
    file_4_occ{counter,1} = set_num;
    file_4_occ(counter,2) = {[num2str(set_num),'_corrprot']};
    file_4_occ{counter,3} = match_prot_log;
    file_4_occ(counter,4) = {s_prot};
    file_4_occ(counter,5) = set_info{i,'GeneSymbol'};
     file_4_occ(counter,6) = set_info(i,5);
      file_4_occ(counter,9) = set_info(i,'PP_Reference');
         
    if (match_prot_log)
        data_store = [data_store; matched_prot_data];

    else
        data_store = [data_store; NaN([1,data_col_num])];
    end
    
    counter = counter+1;
end



end



