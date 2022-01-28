%% organize the phosphodata sets - combine info about where the parent peptide was found

function [sets_info,Mod_Sets_T,Occ_Data] = ...
    PreOcc_Major_Two_part1(filename_Mod,filename_Parent,...
    Channel_Names,ModNonDataVars)

mod_T = readtable(filename_Mod,'format','auto'); % this should be the baciq data
mod_data = mod_T(:,Channel_Names);

phos_sets = mod_T{:,1};
%phos_forms_dat = mod_T{:,2};

[u_p_sets,~,u_p_ic] = unique(phos_sets,'stable'); % set #s are in the second column
sets_info = zeros(length(phos_sets),2);

ModSets_info = [];
ModSets_data = [];
Set = [];
Occ_Data = [];

parent_T = readtable(filename_Parent,'format','auto');
parent_sets_str = string(parent_T{:,'Set'});

parent_sets = zeros(length(parent_sets_str),1);
for i = 1:length(parent_sets_str)
    spl_name = split(parent_sets_str(i),'_');
    parent_sets(i) = str2double(spl_name{end});
end

parent_data = parent_T(:,Channel_Names);
OCC_FLAG = [];

cnt = 0;

for i = 1:length(u_p_sets) % for each set
    
    ModSets_info = [ModSets_info; mod_T(u_p_ic==i,ModNonDataVars)];
    ModSets_data = [ModSets_data; mod_data(u_p_ic==i,:)];
    Set = [Set; u_p_sets(i)*ones(sum(u_p_ic==i),1)];
    
    cnt = cnt+1;
    
    ind_parent = find((parent_sets) == u_p_sets(i));
        
        if (ind_parent)
            
            sets_info(cnt,1) = 1;
            
            Occ_Data = [Occ_Data; parent_data{ind_parent,:}];
            OCC_FLAG = [OCC_FLAG ; ones(sum(u_p_ic==i),1)];
        else
            OCC_FLAG  = [OCC_FLAG ; zeros(sum(u_p_ic==i),1)];
            
        end
        
    
    cnt = cnt+sum(u_p_ic==i)-1;
    
end

Mod_Sets_T = [array2table(Set) array2table(OCC_FLAG) ModSets_info ModSets_data];

end