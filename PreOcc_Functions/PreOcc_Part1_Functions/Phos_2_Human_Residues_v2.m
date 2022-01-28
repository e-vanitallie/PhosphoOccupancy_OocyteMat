function [table_out,count_all_REFS,count_all_RES] = Phos_2_Human_Residues_v2(mod_peps_info)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

prot_refs = mod_peps_info{:,'Protein_Reference'};
[u_prot_refs,ia_pr,ic_pr] = unique(prot_refs,'stable');

count_all_REFS = length(u_prot_refs);

lin_2_inds = 1:length(u_prot_refs);

phos_positions = mod_peps_info{:,'Phosphosite_Position'};

store_phos_pos = cell(length(u_prot_refs),1);

count_all_RES = 0;

phos_motifs = mod_peps_info{:,'Phosphosite_Motif'};
store_phos_motifs = cell(length(u_prot_refs),1);

for i = 1:length(u_prot_refs)
    
    inds_to_positions = find(ic_pr == i);
    phos_store = [];
    motif_store = {};
    
    for j = 1:length(inds_to_positions)
        
        next_to_add = phos_positions{inds_to_positions(j)};
        split_out = split(next_to_add,';');
        
        next_to_add_M = phos_motifs{inds_to_positions(j)};
        split_out_M = split(next_to_add_M,';');
        
        for k = 1:length(split_out)
            new_num = str2double(split_out{k});
            check_already = ismember(new_num,phos_store);
            
            if (check_already) % if it already there do not add it
            else
                phos_store = [phos_store new_num];
                motif_store = [motif_store split_out_M{k}];
            end
            
        end
        
    end
    
    % once we have gone through all the numbers 
    if (phos_store)
        phos_store = sort(phos_store);
        
        m = 0;
        out_str = '';
        out_str_M = '';
        
        count_all_RES = count_all_RES + length(phos_store);
        
        while m < length(phos_store)
            m = m+1;
            out_str=strcat(out_str,';',num2str(phos_store(m)));
            out_str_M = strcat(out_str_M,';',motif_store{m});
        end
        store_phos_pos{i} = out_str;
        store_phos_motifs{i} = out_str_M;
    else
        store_phos_pos{i} = '';
        store_phos_motifs{i} = '';
    end
end

% but now need to remove the references that are to false discoveries or
% contaminants 
refs_contams = contains(u_prot_refs,'contam');
refs_not_contams = ~refs_contams;

refs_fdr = strncmp(u_prot_refs,'##',2);
refs_not_fdr = ~refs_fdr;

u_refs_keep = all([refs_not_contams refs_not_fdr],2);
keep_inds = lin_2_inds(u_refs_keep);

gene_symbs_F = mod_peps_info(ia_pr(keep_inds),'GeneSymbol');
prot_refs_UF = mod_peps_info(ia_pr(keep_inds),'Protein_Reference');

table_out = [gene_symbs_F prot_refs_UF ...
    cell2table([store_phos_pos(keep_inds) store_phos_motifs(keep_inds)],...
    'VariableNames',{'Phosphosite_Position_s','Motifs'})];



end

