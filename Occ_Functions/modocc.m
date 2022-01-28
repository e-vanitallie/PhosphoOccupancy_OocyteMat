function modocc(data_file,target_ci,bootci_iterations,norm2protein,savefigure,occ_plot_F,directoryname,cluster_toggle,occ_x_label,occ_x_ticks)
% modocc('filename', ALPHA, Iterations, Norm2Protein, SavePlot, ComputeParallel, DataVars)
% 
% An algorithm to calculate the site stoichiometry of a phospho-sites (or other modifications) 
% measured under multiple conditions with quantitative mass spectrometry and provide confidence intervals. 
% The approach is compatible with single and multiply modified peptides measured over an arbitrary number of conditions.
%
% Seven required arguments, which are filename (string), target alpha (between zero and one), 
% number of bootstrap interations, and then logicals input for whether to normalize to the protein value,  
% plotting the data and parallel computing toggle. Finally, a cell
% array of strings specifying the data variables of the input table must be included. 
% 
% To run:
% 
%  1) Download repository files and move to the same filepath of your choice.
% 
%  2) Run script using Matlab, which calls the other necessary function 'call_fit_SVD' included in the download.
% 
%  3) Code input is the relative changes of the unmodified and modified forms of a peptide in a matrix.

% This is the output of the matching function. The general format is:
% 
% [ Unmodified ratios....
%   Form 1 ratios ...
%   .
%   .
%   .
%   Form N ratios ...]
%
%  4) Output is a .csv file of the estimated occupancy, the high and low bounds of the confidence intervals, as well as the original input data. Filename is records the parameters and input file use to generate the data along with the date.  
%  The occupancy value returned in the table with a "_occ" suffix added.
%  The confidence bounds have "_hi" and "_lo" suffixes added. 

%  # Visualizations
% 
%  -The code will output two visualizations in a separate folder:
% 
%  1) Subplot showing the unmodified and modified plot over the conditions and the output occupancy values with shaded confidence intervals in each form.
%     In high quality data, the trends should be roughly reciprocal. The values at each of these species at each conditions are used to perform the regression. 
%     The plots Occupancy trend over conditions for each form with confidence intervals. 
%  
%  2) Fitting plots that the solutions are obtained from. Shows the coordinates defined by the measured species
%  and the fitting lines based on the singular value decomposition. The
%  solution spanned by the last principle component (dotted line).
% 
% Tips:
% 
% - Requires Matlab 2014 or later. Pseudocode to assist in nonmatlab implementations is included in Presler et al, 2017.
% - For initial runs, limit the amount of bootstrapping iterations to 10-100 in Section 2, as the code is can be quite slow on most computers for larger datasets when using the recommended 1,000-10,000 iterations.

%% 1) Data import 

%parameter flags
%if ischar(phospho_file)==0
%    error('phospho_file parameter must be characters')
%end

if isnumeric(bootci_iterations)==0
    error('bootci_iterations parameter must be numerical')
end 

if target_ci>1 || target_ci<0
    error('target_ci parameter must be between 0 and 1')
end

if savefigure~=1 && savefigure~=0
    error('plot_data must be logical (0 or 1)')
end

%if cluster_toggle~=0 && cluster_toggle~=1
%    error('cluster_toggle parameter must be logical (0 or 1)')
%end

FullData = readtable(data_file, 'ReadVariableNames', true, 'Delimiter', ',','format','auto');

% for phospho file
%ModNamesRaw = importdata(ModColumnFile);
DataVars = {};
LabelVars = {};
ImportDataVarTypes = {};
for i = 1:size(FullData.Properties.VariableNames,2)
    if regexp(FullData.Properties.VariableNames{i}, 'Chan', 'once') == 1
        DataVars = [DataVars, FullData.Properties.VariableNames{i}];
        ImportDataVarTypes = [ImportDataVarTypes, 'double'];
    else
        LabelVars = [LabelVars, FullData.Properties.VariableNames{i}];
        if strcmp(FullData.Properties.VariableNames{i}, 'set_id')
            ImportDataVarTypes = [ImportDataVarTypes, 'double'];
        else
            ImportDataVarTypes = [ImportDataVarTypes, 'cellstr'];
        end
    end
end

% don't continue if less than two time points found
if size(DataVars, 2) < 2
    error('Please provide more than one time point worth of measurements.')
end

% parse comma separated string to cell array for xlabel ticks
xticks = strsplit(occ_x_ticks, ',');

%DataVars2 = matlab.lang.makeValidName(split(ModNamesRaw{:}, ',')');

%LabelVars2 = {'set_id','site_id','match_protein_log','ProteinID','GeneSymbol','residue_positions'};

% Pulled up the output directory name check from lower in the code to up here
% Don't want to run a bunch of computation and THEN check that the output
% directory already exists.
% This method will  make the script die quickly if you already have a
% directory with the specified name

%Set date for saving files
date_format = 'yyyymmdd';
date_label = datestr(now,date_format);

%if savefigure

    directory = pwd;
    newdirectory  = [directory,filesep, directoryname, '_Output_',date_label];

    if exist(newdirectory, 'dir')
        error('Please choose another output directory name. The specified directory name already exists.')
    end 

    mkdir(newdirectory)
    print_path=[newdirectory,filesep];
    
    if savefigure
    
        newdirectory_fits = [newdirectory,filesep,'Fit_plots',filesep];
        mkdir(newdirectory_fits)
    end
    
    if occ_plot_F
        newdirectory_occ = [newdirectory,filesep,'Occupancy_plots',filesep];
        mkdir(newdirectory_occ)
    end
    
%end

%Automated labeling of filename with bootstrapping parameters
bootstrap_num_label = num2str(bootci_iterations);

if norm2protein
    protnorm_label = 'ProtNorm_';
else
    protnorm_label = '';
end

% [filepath,name,ext] = fileparts(filename) 
% working off of phospho_filename here
% use this for constructing output filename
%[~,filename_basename,~] = fileparts(phospho_file);
[~, filename_basename,~] = fileparts(data_file);
stoich_table_filename = [newdirectory,filesep, date_label,'_Occupancy_trends_',protnorm_label,'ci_',num2str(target_ci*100),'_',bootstrap_num_label,'x_from_',filename_basename,'.csv'];
exit_codes_table_filename = [newdirectory,filesep, date_label,'_Occupancy_trends_',protnorm_label,'ci_',num2str(target_ci*100),'_',bootstrap_num_label,'x_from_',filename_basename,'exit_codes.csv'];



% iterate on the list of unique sets
tableOfTables = {};

exitCodesTableOfTables = {}; 

% create the input structure (one set at a time per loop iteration)
sets = unique(FullData(:,"set_id"));

if cluster_toggle>128
    c = parcluster;
    c.AdditionalProperties.WallTime = '0-8';
    c.AdditionalProperties.QueueName = 'short';
    c.AdditionalProperties.MemUsage = '1000';
elseif cluster_toggle~=0
    c = parcluster('local');
else
    % if cluster_toggle is not 0 or 1 
    % we don't parcluster
end

% No parallelization 0 
% Parcluster local (your computer, every core) 1
% Parcluster local (specified number of cores) > 1 <128
% Cluster (for O2) - not between 0 and 128, >128


for i = 1:size(sets,1)
    ImportData = FullData(FullData.set_id == sets(i,1).set_id,:);
    % cluster parallelization
    if cluster_toggle > 128
        slurm_job_name_parameter = sprintf('-J set_%d_%s', sets(i,1).set_id, filename_basename);
        c.AdditionalProperties.AdditionalSubmitArgs = slurm_job_name_parameter; 
    % checking if local with specified # of cores
    elseif cluster_toggle > 1 && cluster_toggle <= 128
        c.NumWorkers = cluster_toggle;
    end
    % cluster + local parallelization
    if cluster_toggle>=1 
        % the second paremeter is number of output functions
        sprintf('Submitting job for set %d', sets(i,1).set_id);
        job(i) = c.batch(@runOcc, 5, {ImportData, DataVars, LabelVars, print_path, norm2protein, savefigure,occ_plot_F, bootci_iterations, target_ci, occ_x_label, xticks});
    % no parallelization
    else
        [ExportTable, was_bootci_called, exit_codes, Number_of_Forms, setID] = runOcc(ImportData, DataVars, LabelVars, print_path, norm2protein, savefigure, occ_plot_F, bootci_iterations, target_ci, occ_x_label, xticks);
        tableOfTables{i} = ExportTable;
        
        exit_code_vals = exit_codes{1,:};
        exit_code_var_names = exit_codes.Properties.VariableNames;
      %  exitCodesTableOfTables{i} = array2table([(setID), Number_of_Forms,  was_bootci_called,  join(string(table2cell(exit_codes)), ",")], 'VariableNames', {'set_id', 'num_phospho_forms', 'boot_ci_called', 'Exit_codes_(underdetermined,#_uncorrected,form_data_identical)'});
        exitCodesTableOfTables{i} = array2table([setID Number_of_Forms was_bootci_called exit_code_vals], 'VariableNames', [{'set_id', 'num_phospho_forms', 'boot_ci_called'},exit_code_var_names]);
   
    end
end


if cluster_toggle >=1
    for i = 1:size(sets,1)
       % job(i).wait
       wait(job(i),'finished')
    end
    for i = 1:size(sets,1)
        if job(i).fetchOutputs{2} == 1
            tableOfTables{i} = job(i).fetchOutputs{1};
        end
        % "Set #" "# phospho forms" "Boot ci called?" "Flags"
        exit_code_T = job(i).fetchOutputs{3};
        exit_code_vals = exit_code_T{1,:};
        exit_code_var_names = exit_code_T.Properties.VariableNames;
        exitCodesTableOfTables{i} = array2table([job(i).fetchOutputs{5} job(i).fetchOutputs{4} job(i).fetchOutputs{2} exit_code_vals], 'VariableNames', [{'set_id', 'num_phospho_forms', 'boot_ci_called'},exit_code_var_names]);
    end
end

fullTable = vertcat(tableOfTables{:});

fullExitCodesTable = vertcat(exitCodesTableOfTables{:});



%% 6) Consolidate and export data (not fully functioning at moment) 
% 
% Returns:      
%   The code exports a table of  occupancies with the upper and lower bound 
% of the confidence intervals at a given cutoff, if specified. 

%save table
writetable(fullTable,stoich_table_filename);

writetable(fullExitCodesTable, exit_codes_table_filename, 'Delimiter', ',');

end
