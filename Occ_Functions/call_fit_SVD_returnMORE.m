function [out_info, Smatrix, Vmatrix, normal_vector] = ...
    call_fit_SVD_returnMORE(input_matrix,MC_flag)



if (MC_flag == 0)
     MeanCenteredData = input_matrix - repmat(mean(input_matrix,1),length(input_matrix),1);    
        u_check = sum(sum(logical(MeanCenteredData)));
        
        if (u_check == 0)
            occupancy_trend = 100*round(rand(size(MeanCenteredData,2),1),0);
            imposible_fit = 2;
             out_info = [imposible_fit; occupancy_trend];
            return
        end
        
        input_matrix = MeanCenteredData;        
end

    
[~, Smatrix, Vmatrix] = svd(input_matrix);

normal_vector = Vmatrix(:,end);

%Calculates occupany from the orthogonal vector (i.e., normal) to the plane
%defined by the experimental values. We find the vector whose elements sum
%to 1 (similar to a stochastic or probability vector) to satisfy the conservation law. 
%We multiply by 100 to convert to percent occupancy. 

nv_r3 = round(normal_vector,3); % ROUND values between (-0.0005,0) to 0

occupancy_trend_S = normal_vector./repmat(norm(normal_vector,1),length(normal_vector),1)*100;
occ_T_r3 = nv_r3./repmat(norm(nv_r3,1),length(nv_r3),1)*100;

if ( norm(occupancy_trend_S,1) == abs(sum(occupancy_trend_S)) ) 
    % these are only equal if all values are either positive or negative
    occupancy_trend = abs(occupancy_trend_S);
    imposible_fit = 0;
    out_info = [imposible_fit; occupancy_trend];
    return

elseif ( norm(occ_T_r3,1) == abs(sum(occ_T_r3)) )
    % check if the rounding alows us to find the right answer
    
    if ( norm(occ_T_r3,1) ~= 100 )
        % scale occ_T_r3
        occ_T_S = occ_T_r3*(100/norm(occ_T_r3,1));
        occupancy_trend = abs(occ_T_S);
        imposible_fit = 1;
    else
        occupancy_trend = abs(occ_T_r3);
        imposible_fit = 1;
    end
    
    out_info = [imposible_fit; occupancy_trend];
    return
    
else % if neither of the two other conditions are met, return 0s and 100s
    
    occupancy_trend = 100*round(rand(size(normal_vector,1),size(normal_vector,2)),0);
    impossible_fit = 2;
    out_info = [impossible_fit; occupancy_trend];
    
end
     
   
end

