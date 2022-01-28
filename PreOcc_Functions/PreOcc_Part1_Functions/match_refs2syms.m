function out_1=match_refs2syms(Refs,db1)

% Refs is the cell array of the references that we are going to match to
% gene symbols

% db1 is the list of all the gene symbols we would like to match to 

% out_1 is vector of the same length as "Refs" where there are 0s if the
% reference was not found in db1 and there is the value of the reference in
% db1 which allows for finding corresponding gene symbols in the next step

out_1 = zeros(1,length(Refs));

    for i = 1:length(Refs)
    
    s = Refs{i};
    l = length(s);
    f = strncmp(s,db1,l);
    
    ind = find(f);
    if (isempty(ind))
    elseif length(ind)>1
        for j = 1:length(ind)
            f_rev = strncmp(db1{j},s,length(db1{j}));
            if (f_rev)
                out_1(i) = ind(j);
            end
        end
    else
        out_1(i) = ind;
    end
    
    end
end