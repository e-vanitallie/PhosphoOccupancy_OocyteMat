function [exit_codes, was_bootci_called, ExportTable, setID] = ...
    computeOcc(ImportData, norm2protein, DataVars, print_path, savefigure, ...
    bootci_iterations, target_ci, LabelVars,...
    occ_x_label,occ_x_ticks)

FP_color_MATRIX = [1.0000    1.0000    1.0000;...
    1.0000         0         0; ...
    0.9290    0.6940    0.1250; ...
    1.0000    1.0000         0; ...
         0    1.0000         0; ...
         0    1.0000    1.0000; ...
         0    0.4470    0.7410; ...
         0         0    1.0000; ...
    1.0000         0    1.0000; ...
    0.7500    0.7500    0.7500; ...
         0         0         0];

low_p_thresh = 0.02;

% computeOcc.m does occupancy calculation, bootstrapping, normalization, fit plotting

% we define a set of exit codes to denote status for the given set
% being processed. We return the corresponding exit code depending on
% which case is handled.
%
% Exit codes (in the order they are implemented in the function):
% 0 - we get to the end with no issues
% 1 -

%Sets consistent random number set
rng('default')

setID = ImportData.set_id(1,1);
Has_Protein = ImportData.Match_Protein;

%No protein index, 0 if corrprot in site_ids, otherwise 1
% corrprot is data for normalization

FullsetData = ImportData{:,DataVars};

if isnumeric(FullsetData)==0
    error('Input matrix must be numerical')
end % end of isnumeric if block

% everything except for the norm protein line (why include parent?)
ModsetData = FullsetData(1:(end-1),:);
Orig_ModsetData = ModsetData; % save in case need for plotting removed points

% just the norm protein
ProtData = FullsetData(end,:);

[NumForms, NumConditions] = size(ModsetData);

%Analysis cannot proceed if the system is underdetermined.
underdetermined_flag = 0;
if NumForms > NumConditions
    
    underdetermined_flag = 1;
    warning('System underdetermined- the number of modified forms exceeds the number of measured conditions in input matrix')
    %TODO: return some sort of indicative exit code
end  % end of underdetermined (NumForms > NumConditions) if block

%Used to make fitting figures if set to 1 above
if savefigure
    fitting_plot = figure('visible','off');
    
    dim_sp_1 = round(sqrt(length(Orig_ModsetData)));
    dim_sp_2 = round(length(Orig_ModsetData)/dim_sp_1);
    if (dim_sp_1*dim_sp_2)<length(ModsetData)
        dim_sp_2 = dim_sp_2 + 1;
    end % end of dim_subplot_second if
    %fitting_plots = figure('visible','on');
end % end of save figure

% for keeping track of the number that have not been corrected
% initialize to 0 if the main analysis section is not entered
total_num_uncorrected = 0;
% flag to retain whether or not we called bootci
was_bootci_called = 0;

%%%%%%  Enter main analysis section (if critera from above are met): %%%%
if (underdetermined_flag == 0)
    
    total_num_uncorrected = NumConditions;
    NumConditions_orig = NumConditions;
    conditions_eval = 1:1:NumConditions;
    should_plot = 0;
    
    %Code will fail if there are zeros in channels, since it may divide by zero later on.
    %Replace zeros with 1E-9 for all forms
    ModsetData(ModsetData==0) = 1E-9;
    ProtData(ProtData==0) = 1E-9;
    ProtDataOrig = ProtData;
    modset_NOdetrend = ModsetData;
    
    if norm2protein == 1 && Has_Protein(1) == 1
        
        low_prot_check = logical(ProtData <= low_p_thresh);
        %save_rel = ModsetData; % save the original data before deterending
        %save_prot = ProtData; % save the original protein data
        
        if (sum(low_prot_check)) % if there are values w/ too low protein
            
            total_num_uncorrected = sum(~low_prot_check); % have to re-set total num uncorrected
            
            if (total_num_uncorrected < NumForms)
                exit_codes = array2table([1, 0, 0],...
    'VariableNames', {'underdetermined_flag','total_num_uncorrected', 'extra'});
                ExportTable = [];
                return;
            end
            
            ProtData = ProtData(:,~low_prot_check);
            ModsetData = ModsetData(:,~low_prot_check);
            NumConditions = NumConditions - sum(low_prot_check);
            
            conditions_eval = conditions_eval(~low_prot_check);
            
            occ_x_ticks = occ_x_ticks(~low_prot_check);
            FP_color_MATRIX = FP_color_MATRIX(~low_prot_check,:);
            
        end
        
        
        prot_detrend = ProtData; %./repmat(min(ProtData),1,size(ProtData,2));
        % scale where 1 is lowest so that scale is clear for
        % de-bugging but could not use
      
        modset_detrend = ModsetData./repmat(prot_detrend,size(ModsetData,1),1);
        ModsetData = modset_detrend;
        
        
    else
        low_prot_check = ~true(size(ProtData,1),size(ProtData,2));
        %save_rel = ModsetData;
        
    end % end of norm2protein and has_protein check
    
    %preallocate
    OccupancyAllConditions = zeros(NumForms,NumConditions);
    InputMatrixAllConditions = zeros(NumConditions,NumForms,NumConditions);
   
    if (NumForms > 3)
        %disp('hi')
    end
    
    %Calculate optimal solution for time series
    for condition_counter1 = 1:NumConditions
        
        %Normalize to reference timepoint and reshape matrix to NumConditions,NumForms for compatablility with bootci input format
        scaleMatrix = [ModsetData./repmat(ModsetData(:,condition_counter1),1,NumConditions)]';
        % Mean Center Input Matrix
        MeanCenteredData = scaleMatrix - repmat(mean(scaleMatrix,1),length(scaleMatrix),1);
       
        
        u_check = sum(sum(logical(MeanCenteredData)));
        
        if (u_check == 0)
            OccupancyEstimate = 100*round(rand(size(MeanCenteredData,1),1),0);
        else
            
            InputMatrix = MeanCenteredData;
            %Stores input matrix for use in bootstrapping so it's not necessary to repeat code below
            InputMatrixAllConditions(:,1:NumForms,condition_counter1) = scaleMatrix;
           
            %Perform fitting to calculates percent occupancy using the relative data of each form as input
            %Exports the occupancy estimate for each times point (corrected
            %and uncorrects for "impossible" values above or below zero),as
            %well as standard SVD output matrices.
            [OccupancyEstimate, Smatrix, Vmatrix, normal_vector] = ...
                call_fit_SVD_returnMORE(InputMatrix,1);
            %            [OccupancyEstimate, OccupancyEstimateUncorrected, MeanCenteredData, Smatrix, Vmatrix, normal_vector, isCorrected,Umatrix]...
            %                = call_fit_SVD(InputMatrix, FLAG_repeated_points_bstrp);

        end
        
        % decrement total_num_uncorrected.
        % if it has been corrected, we subtract one!
        % if it hasn't been corrected, 0 will be subtracted.
        impossible_check = logical(OccupancyEstimate(1)==2);
        total_num_uncorrected = total_num_uncorrected - impossible_check;
        
        OccupancyAllConditions(:,condition_counter1) = OccupancyEstimate(2:end);
        
        %Visualizing SVD fitting data
        
        % flag to turn off plotting for whatever reason (e.g. if
        % NumForms >3), initially this is set to "yes, we should
        % plot" or 1.

        
        if ( savefigure && logical(NumForms <= 3))
            
            should_plot = 1;
            subplot_title = [occ_x_label,'',occ_x_ticks{condition_counter1}];
            sp = conditions_eval(condition_counter1);
            
            [fitting_plot,p_normal_v,data_p] = plot_SVD_out(fitting_plot,NumForms,[dim_sp_1 dim_sp_2],...
                sp,subplot_title,FP_color_MATRIX,NumConditions,...
                InputMatrix,Vmatrix,Smatrix,normal_vector,round(OccupancyEstimate(2:end),1),impossible_check);
        
       
        end   % end of savefigure if
        
        
        
    end % end of "calculate optimal solution for time series" for loop
%     
    if ( savefigure && should_plot)
        sgtitle(['Occupancy Estimate for each ',occ_x_label,', using complete data'],'FontSize',12)
        fitting_plot.Visible = 'off';
        %figure(fitting_plot)
        % figure = fitting_plot('visible','off');
        if NumForms == 2
            h1 = legend([p_normal_v data_p],[{'Normal Vector'},occ_x_ticks]);
        else
            h1 = legend([p_normal_v data_p],[{'Normal Vector','2D Plane'},occ_x_ticks]);
        end
        p = get(h1,'Position');
        p(1) = 0.8;
        p(2) = 0.08;
        set(h1,'Position',p)
    end % end of if for saving figures
    
    if savefigure && (total_num_uncorrected > 0)
        print(fitting_plot,'-dpng',[print_path,'Fit_plots/','OccEst_Set_',num2str(setID(1)),'_fits'])
        close(fitting_plot)
    elseif should_plot == 1 && savefigure
        print(fitting_plot,'-dpng',[print_path,'Fit_plots/','NoOcc_Set_',num2str(setID(1)),'_fits'])
        close(fitting_plot)
    elseif savefigure
        close(fitting_plot)
    end% end of if for saving figures
    
    
    % 4) Determine CIs with boostrapping
    
    
    if (total_num_uncorrected > 0)
        
        was_bootci_called = 1;
        %preallocate
        Fit_confidence_lower_SVD = zeros(NumForms, NumConditions);
        Fit_confidence_higher_SVD = zeros(NumForms, NumConditions);
       
        %protoype for stepping through particular sets...
        for ConditionsCounter_CIs = 1:NumConditions
            
            input_matrix_CI = InputMatrixAllConditions(:,:,ConditionsCounter_CIs);
            MC_flag = 0;
            
        % Use BOOTCI and then you can get bias corrected confidence
        % intervals 
            ci_occupancies_SVD = bootci(bootci_iterations,{@call_fit_SVD,input_matrix_CI,MC_flag},...
                'alpha',(1-target_ci),'type','bca','options',...
                statset('UseParallel',0,'UseSubstreams',0));
           ci_occupancies_SVD = ci_occupancies_SVD';
 Fit_confidence_lower_SVD(:,ConditionsCounter_CIs) = ci_occupancies_SVD(:,1);
 Fit_confidence_higher_SVD(:,ConditionsCounter_CIs) = ci_occupancies_SVD(:,2);

        % Use BOOTSTRAP and you can return the actually values of the
        % distribution 
        
                        
        end % end of "stepping through particular sets" section
        
    end
end

%     Store all the occupancies and their confidence bounds for each site
%     We don't need to track positions using insert_positions, as we're
%     creating independent data structures (e.g. Optimal_occupancy_stored)
%     with each time computeOcc is called
%     Previously, we'd be appending to the same data structures for all
%     sets and would need to figure out the offset, essentially.
%     insert_positions = stored_data_counter:((stored_data_counter-1)+size(OccupancyAllConditions,1));


if (was_bootci_called)
    
    % NEED TO stich back matrix from removing low protein points
    Optimal_occupancy_stored = NaN(size(ModsetData,1),NumConditions_orig);
    CI_lower_stored = NaN(size(ModsetData,1),NumConditions_orig);
    CI_higher_stored = NaN(size(ModsetData,1),NumConditions_orig);
    relative_trends_stored = NaN(size(ModsetData,1),NumConditions_orig);
    modset_NOdetrend_stored = NaN(size(ModsetData,1),NumConditions_orig);
    
    Optimal_occupancy_stored(:,~low_prot_check) = OccupancyAllConditions;
    CI_lower_stored(:,~low_prot_check) = Fit_confidence_lower_SVD;
    CI_higher_stored(:,~low_prot_check) = Fit_confidence_higher_SVD;
    
    relative_trends_stored(:,~low_prot_check) = ModsetData; % save the detrended results !!
    modset_NOdetrend_stored(:,:) = modset_NOdetrend;
    prot_data_out = repmat(ProtDataOrig,size(ModsetData,1),1);
    
    Combined_data = [Optimal_occupancy_stored, CI_higher_stored, CI_lower_stored, relative_trends_stored,...
        modset_NOdetrend_stored,prot_data_out];
    
    % labels_opt contains values such as:
    labels_opt = cell(1,length(DataVars));
    labels_hi = cell(1,length(DataVars));
    labels_lo = cell(1,length(DataVars));
    labels_orig = cell(1,length(DataVars));
    labels_prot = cell(1,length(DataVars));
    
    for RenameCounter = 1:length(DataVars)
        labels_opt(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_occ']);
        labels_hi(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_hi']);
        labels_lo(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_lo']);
        labels_orig(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_orig']);
        labels_prot(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_prot']);
    end
    
    % combine all sets of labels
    Data_labels_combined = [labels_opt,labels_hi,labels_lo, DataVars,labels_orig,labels_prot];
    
    % put the occupancy, hi, low, and ModsetData (parent + phospho, NO
    % NORM protein line) into CombinedDataTable
    
    CombinedDataTable = array2table(round(Combined_data, 5),'VariableNames',Data_labels_combined);
    labeltable = ImportData(1:NumForms,LabelVars);
    
    ExportTable = [labeltable,CombinedDataTable];
    
else
    ExportTable = [];
end

% return the exit codes, yipee
extra_flag = 0;
exit_codes = array2table([underdetermined_flag, total_num_uncorrected, extra_flag],...
    'VariableNames', {'underdetermined_flag','total_num_uncorrected', 'extra'});

% function end
end