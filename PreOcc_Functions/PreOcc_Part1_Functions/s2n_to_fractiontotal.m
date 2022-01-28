function table_out = s2n_to_fractiontotal(file_in,Chan_In)
data_A = readtable(file_in,'format','auto');
data_N = data_A{:,Chan_In};

sum_all = sum(data_N,2);
data_frac = round(data_N./repmat(sum_all,1,length(Chan_In)),3);

table_out = [data_A(:,1) array2table(data_frac,'variablenames',Chan_In)];

end