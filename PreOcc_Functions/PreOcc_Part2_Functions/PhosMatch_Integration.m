function comb_phos_table = PhosMatch_Integration(info_matching, phos_MATCH_filename)

human_match_out = readtable(phos_MATCH_filename,'format','auto');
xen_refs_matched = human_match_out{:,1};
xen_residue_matched = human_match_out{:,4};

store_match = zeros(size(info_matching,1),1);

for i = 1:height(info_matching)
    
    xen_ref_M = info_matching{i,5};
    xen_ref_res = split(info_matching{i,6},';');
    
    check_ref = find(strcmp(xen_ref_M, xen_refs_matched));
    
    if sum(check_ref)
        
        for j = 1:length(xen_ref_res)
            
            [res_check,res_loc] = ismember(str2double(xen_ref_res{j}),xen_residue_matched(check_ref));
            
            if (res_check)
                
                store_match(check_ref(res_loc)) = i;     
                
            end
        end
        
    end 
end

comb_phos_table = table('size',[height(info_matching),6],...
    'variabletypes',["cellstr","cellstr","cellstr","cellstr","cellstr",...
    "cellstr"]);
comb_phos_table.Properties.VariableNames = ...
    {'Human_Reference','Match_Code','Human_Residue','Human_Motif','Motif_Score',...
'Human_LT_Info'};

store_match = store_match(logical(store_match));

for i = 1:length(store_match)
    
    comb_phos_table{store_match(i),1} = ...
        human_match_out{i,comb_phos_table.Properties.VariableNames{1}};
    
    if ~ismissing(comb_phos_table{store_match(i),2})
        
        comb_phos_table{store_match(i),2} = ...
            cellstr(strcat(comb_phos_table{store_match(i),2},";",...
            num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{2}})));
        
         comb_phos_table{store_match(i),3} = ...
            cellstr(strcat(comb_phos_table{store_match(i),3},";",...
            num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{3}})));
        
          comb_phos_table{store_match(i),4} = ...
            cellstr(strcat(comb_phos_table{store_match(i),4},";",...
            (human_match_out{i,comb_phos_table.Properties.VariableNames{4}})));
        
          
         comb_phos_table{store_match(i),5} = ...
            cellstr(strcat(comb_phos_table{store_match(i),5},";",...
            num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{5}})));
        
         comb_phos_table{store_match(i),6} = ...
            cellstr(strcat(comb_phos_table{store_match(i),6},";",...
            num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{6}})));
        
%          comb_phos_table{store_match(i),7} = ...
%             cellstr(strcat(comb_phos_table{store_match(i),7},";",...
%             num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{7}})));
%         
        
        
    
    else
        comb_phos_table(store_match(i),2) = ...
            {num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{2}})};
        
          comb_phos_table(store_match(i),3) = ...
            {num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{3}})};
  
        comb_phos_table{store_match(i),4} = ...
            human_match_out{i,comb_phos_table.Properties.VariableNames{4}};
  
        
       comb_phos_table(store_match(i),5) = ...
            {num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{5}})};
  
        
        comb_phos_table(store_match(i),6) = ...
            {num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{6}})};
  
        
%         comb_phos_table(store_match(i),7) = ...
%             {num2str(human_match_out{i,comb_phos_table.Properties.VariableNames{7}})};
%   
%         
        
    end
    
end

comb_phos_table.Properties.VariableNames = ...
    {'PP_Reference','PP_Match','Human_Residue','Human_Motif','Motif_Score',...
'Human_LT_Info'};



end
