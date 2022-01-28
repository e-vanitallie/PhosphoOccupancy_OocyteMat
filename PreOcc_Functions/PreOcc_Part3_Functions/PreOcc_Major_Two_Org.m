function [MatchedPhos_T] = ...
    PreOcc_Major_Two_Org(col_nums_out,Channel_Names,filename_Mod,...
    ModNonData,filename_Parent,...
    prot_proteinIDs,prot_data,mod_occ_out_FixedNames,...
    occ_table_var_types_extra)


ModNonDataVars = ModNonData;

[sets_indexing,Mod_Sets_T,Parent_Data] = ...
PreOcc_Major_Two_part1(filename_Mod,filename_Parent,...
Channel_Names,ModNonDataVars);

set_info = Mod_Sets_T(logical(sets_indexing(:,1)),...
    ModNonDataVars);

phos_info = Mod_Sets_T(logical(Mod_Sets_T.OCC_FLAG),...
    [{'Set','P_Form','Phosphosite_Position','Phosphosite_Motif',...
    'PP_Reference','Human_Residue','Human_Motif','MoreInfoFlag'}]);

occ_phos_data = Mod_Sets_T{logical(Mod_Sets_T.OCC_FLAG),Channel_Names(1:col_nums_out)};

[file_4_occ,data_store] = ...
    PreOcc_Major_Two_part2(set_info,phos_info,occ_phos_data,...
    Parent_Data,prot_proteinIDs,prot_data,occ_table_var_types_extra);

% OUTPUTS from make_MatchedSets 
% (all of these outputs are the same height so that they can be put togehter into output files) 
% 
% 
% file_4_occ -- has the non-DATA information for the mod_occ file
%
% data_store -- has the DATA for the mod_occ file


file_4_occ.Properties.VariableNames = mod_occ_out_FixedNames;

MatchedPhosSets_data = array2table(data_store,'VariableNames',Channel_Names);

MatchedPhos_T = [file_4_occ MatchedPhosSets_data];


end