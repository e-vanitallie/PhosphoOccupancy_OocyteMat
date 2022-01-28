function [out_p1,fraction_detrend] = PreOcc_Major_One_Detrend(madphos_T,phos_info_cols,D_col_nums,...
    Prot_Info_T,Protein_Data)

madphos_T_info = madphos_T(:,phos_info_cols);

Channel_Names_Phos = cell([1,D_col_nums]);
Channel_Names_Protein = cell([1,D_col_nums]);
Channel_Names_Phos_Detrend = cell([1,D_col_nums]);

for i = 1:D_col_nums
    
    Channel_Names_Phos{i} = strcat('Chan_',num2str(i));
    Channel_Names_Protein{i} = strcat('Chan_',num2str(i),'_Prot');
    Channel_Names_Phos_Detrend{i} = strcat('Chan_',num2str(i),'_Detrend');

end

data_to_detrend = madphos_T{:,Channel_Names_Phos};
flag_prot = zeros(size(data_to_detrend,1),1);
indexing_inds = zeros(size(data_to_detrend,1),1);
Prot_Refs4Detrend = Prot_Info_T{:,'PROT_ProteinReference'};
store_detrend = zeros(size(data_to_detrend,1),D_col_nums);

for i = 1:size(data_to_detrend,1)
    find_P = find(strcmp(char(madphos_T_info{i,'Protein_Reference'}),Prot_Refs4Detrend));
    
    if (find_P)
        flag_prot(i) = 1;
        indexing_inds(i) = find_P;
        prot_det = Protein_Data(find_P,:);
        prot_det(~logical(prot_det)) = 1e-6;
        detrend_dat = data_to_detrend(i,:)./prot_det;
        store_detrend(i,:) = detrend_dat./sum(detrend_dat);
    else
        store_detrend(i,:) = NaN;
    end
end
flag_prot = logical(flag_prot);
fraction_detrend = sum(flag_prot)/length(flag_prot);

% organize output table

out_p1 = [madphos_T_info(flag_prot,:) ...
    array2table(round(store_detrend(flag_prot,:),3),'VariableNames',Channel_Names_Phos_Detrend) ...
    array2table(Protein_Data(indexing_inds(flag_prot),:),'VariableNames',Channel_Names_Protein) ...
    madphos_T(flag_prot,Channel_Names_Phos)];
end

