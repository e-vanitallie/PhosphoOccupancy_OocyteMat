function [ExportTable, was_bootci_called, exit_codes, Number_of_Forms, setID] = runOcc(ImportData, DataVars, LabelVars, print_path, norm2protein, savefigure, occ_plot_F, bootci_iterations, target_ci, occ_x_label, occ_x_ticks)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    Number_of_Forms = size(ImportData, 1) - 1;
    
    % call computeOcc.m for occupancy calculation, bootstrapping,
    % normalization, fit plotting, etc.
    
    [exit_codes, was_bootci_called, ExportTable, setID] = ...
        computeOcc(ImportData, norm2protein, DataVars, print_path, ...
        savefigure, bootci_iterations, target_ci, LabelVars,...
        occ_x_label,occ_x_ticks);

%     [exit_codes, was_bootci_called, ExportTable, setID] = ...
%         computeOcc_4(ImportData, norm2protein, DataVars, print_path, ...
%         savefigure, bootci_iterations, target_ci, LabelVars,...
%         occ_x_label,occ_x_ticks);
    
%      [exit_codes, was_bootci_called, ExportTable, setID] = ...
%         computeOcc_EXAMPLE_P(ImportData, norm2protein, DataVars, print_path, ...
%         savefigure, bootci_iterations, target_ci, LabelVars,...
%         occ_x_label,occ_x_ticks);

    % if savefigure is toggled, then we are calling plotocc to save
    % occupancy plots
    
    
    if (was_bootci_called && occ_plot_F) 
        plotOcc_wAll_Info(DataVars, print_path, occ_x_label, occ_x_ticks, ExportTable); % plotting function call
    end
end

