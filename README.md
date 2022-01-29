# PhosphoOccupancy_OocyteMat

This is the READ_ME file for the code for estimating phospho-occupancy from protein and phospho-peptide TMT quantified data processed through Gygi core. 

Also check out “PhosOcc_Schema.png” for a more information about what each script actually does and a visual representation of how things fit together. 

Elizabeth S. Van Itallie; January 28th, 2022

1) Define what you experiment identifier is going to be (i.e. “exp_str”). It is very important to be consistent. For the Oocyte Maturation project, I used “Oocyte_Mat.”

2) Make sure you have the correct folders with the correct inputs in the working directory:

	A. Occ_Functions/

		This folder has the functions that you need to estimate occupancy locally and on the server. They do not change depending on the experiment.

		The functions inside are:
			modocc.m
			runOcc.m
			computeOcc.m
			call_fit_SVD_returnMORE.m
			call_fit_SVD.m
			plot_SVD_out.m
			plotcube.m
			plotOcc_wAll_Info.m 

	B. PreOcc_Functions/ 
		
		This folder has the functions that you need to go from the CORE output files to the input for estimating occupancy. These functions also allow you to generate 1) de-treneded relative phospho-data and 2) protein data with parent peptides removed. There are three subfolders within this folder. 

		PreOcc_Part1_Functions/
		
			These are the functions that are called by “Exp_Str”_PreOcc_Script_Part1.mlx
			
				PreOcc_Major_One_Phos.m
				protein_prepare_p1_v1.m
				phos_prepare_p1_v4.m
				Phos_2_Human_Residue_v2.m
				s2n_to_fractiontotal.m
				match_refs2syms.m

		PreOcc_Part2_Functions/
		
			This is the functions that are called by “Exp_Str”_PreOcc_Script_Part2.mlx

				PhosMatch_Integration.m 
		
		PreOcc_Part3_Functions/

			These are the functions that are called by “Exp_Str”_PreOcc_Script_Part3.mlx
			
				PreOcc_Major_One_Unmod.m
				PreOcc_Major_One_Detrend.m 
				PreOcc_Major_Two_Org.m
				PreOcc_Major_Two_part1.m
				PreOcc_Major_Two_part2.m

	C. “exp_str”_PreOcc_Occ_Scripts/
		
		This folder contains the live MATLAB scripts that call the functions and organize the inputs and outputs so that occupancy can be estimated and the other outputs can be generated. Since these scripts all have variable definition sections at the beginning of the script they will change with each experiment. The scripts should be renamed to start with the “exp_str.” There are three scripts for the “Pre Occupancy” (PreOcc) part of the pipeline and there is one script if you want to estimate “Occupancy” (Occ) locally. 

		The scripts in this folder are: 
		
			“exp_str”_PreOcc_Script_Part1.mlx
			“exp_str”_PreOcc_Script_Part2.mlx
			“exp_str”_PreOcc_Script_Part3.mlx
			“exp_str”_Occ_Script.mlx

	D. “exp_str”_InputFiles_PreOccupancy/
	
		This folder contains the input files that are necessary for the non-human matching part of the occupancy estimation pipeline. 
		
		The files that need to be in this folder are: 
		
			The FASTA file that the mass spectrometer spectra was searched against to identify peptides. 
				For  the oocyte maturation experiment the FASTA file is: “XENLA_9.2_Xenbase_plusMT.pep.fa” 

			The .xlsx file that has the result of reciprocal blast of the Xenopus protein models and human protein models to assign human gene symbols to the Xenopus protein models. 
				For the oocyte maturation experiment the file is: “210113_xenla_v9p2_humangn.xlsx”

			Three different output files that come from Gygi CORE. 
			
			1 - The peptide level protein quantification file. This file should end “protein_quant_with_peptides_NUMBER” where NUMBER is the Protein Quant identification number. 
			
			For oocyte maturation: “OocyteMat_Haas2018ProteinPhos_11plex_siteQuant_5927_peptide_export.csv”

			2 - The peptide level phosphorylation quantification file. This file should end with “siteQuant_NUMBER_peptide_export” where NUMBER is the Site Quant identification number. 

			For oocyte maturation: “OocyteMat_Haas2018ProteinPhos_11plex_siteQuant_5927_peptide_export.csv”

			3 - The composite site phosphorylation quantification file. Composite sites are sites that have multiple phosphorylation residues where at least one of the residues is measured as phosphorylated alone or as a different combination. This file is necessary for accessing the composite site identification strings. We will not use the quantitation data in this file. The file should end with “siteQuant_NUMBER_compositeSite_export” where NUMBER is the Site Quant identification numbers. 
		
	For the oocyte maturation: “OocyteMat_Haas2018ProteinPhos_11plex_siteQuant_5927_compositeSite_export.csv”			
	** The NUMBERS for the two SiteQuant files should be the same!
	** These files are downloaded from CORE as .tsv files. Open them in a text editor and replace the tabs with commas and save as .csv files. 

	E. “exp_str”_Phos_Matching_Python/
		
		This folder contains the python functions, python jupyter notebooks, and input files necessary to identify if the measured Xenopus residues align to homologous human phosphorylation sites and motifs. ** Currently the code is hardcoded MAC operating systems. ** 

		
		The Python3 scripts are: 

			“exp_str”_Matching_Script.py
			step0_match.py
			step1_match.py
			step2_match.py
			step3_match.py
			step4_match.py

		The Python3 Jupyter notebooks are: 

			“exp_str”_Motif_Scorng.ipynb
			“exp_str”_HumanMatches_to_HumanInfo.ipynb 

		The input files are: 
			human-phosphorite-fastas.fasta
			Phosphorylation_site_dataset_032020 
		
		** If you have downloaded this code from GITHUB you need to download your own versions of these files from 
https://www.phosphosite.org/staticDownloads.action  

3) Open MATLAB (this code was written with R2020a) and set your working path to the folder that has all of the folders outlined in “2.” DO NOT add the subfolders to your path. The addition of subfolders is automated. 

4) Now execute the different scripts to generate the input for occupancy  

	A.  First …. OocyteMat_PreOcc_Script_Part1.mlx

		DEFINE THE INPUTS: 
			- “exp_str” — the identifier for the experiment
			- “input_file_folder” - The name of the folder with the input files (e.g. protein references)
			- The names of the input files in the input file folder: 
				“filename_reference” 
				“prot_peps_search_file” 
				“phos_peps_search_file”
				“phos_CompositeInfo_search_file”
			- “delim_core” - the delimiter for the search files 	
			- “protein_measurement_chans” and “phos_measurement_chans” -  The names of the protein and phospho peptide TMT signal channels. The can change depending on how the files are outputted from Core. Make sure that you examine the input files closely to make you designate the right names. The columns should end in “sn.” If the column names in the input files start with a number, then MATLAB will automatically add an “x” in front of the numbers. 
			- “phos_search_class_Col,” “prot_search_class,” “phos_search_class_Col,” “phos_search_class” -
Since the protein and phospho- data are searched on the same protein map to control for the total FDR, the protein and phospho files have quantification information for both the protein and phospho searches. The information from the two different experiments is designated by the “Class” parameter. You need to tell the script what the “Class” columns names are for the protein and phospho files, and what the class identifiers are for the protein and phospho-files.
			- “prot_peps_analyze_file” and “phos_peps_analyze_file”  — The names of the output files that will have the protein and phosphor-peptide information without the irrelevant other “class”
			- “data_col_num” — The number of TMT channels that were used for the experiment. 

		The input for the oocyte maturation project are below. 

	B. Press Run. After is finishes check that a new folder was created: “date_str”_”exp_str”_PreOcc_Files_Part1. There should be 12 output files inside of it.   

% Experiment Identification String
exp_str = 'OocyteMat';

% Input File Folder and input file names 
input_file_folder = ['OocyteMat_InputFiles_PreOccupancy',filesep];

filename_reference = [input_file_folder,'XENLA_9.2_Xenbase_plusMT.pep.fa'];

prot_peps_search_file = [input_file_folder,...
    'OocyteMat_Haas2018ProteinPhos_11plex_protein_quant_with_peptides_26369.csv'];
phos_peps_search_file = [input_file_folder,...
    'OocyteMat_Haas2018ProteinPhos_11plex_siteQuant_5927_peptide_export.csv'];

phos_CompositeInfo_search_file = [input_file_folder,...
    'OocyteMat_Haas2018ProteinPhos_11plex_siteQuant_5927_compositeSite_export.csv'];

delim_core = ',';

% Column names for the TMT quantification for the protein and
% phospho-peptides 
protein_measurement_chans = {'Channel_126_sn','Channel_127n_sn','Channel_127c_sn','Channel_128n_sn',...
    'Channel_128c_sn','Channel_129n_sn','Channel_129c_sn','Channel_130n_sn',...
    'Channel_130c_sn','Channel_131_sn','Channel_131C_sn'};

phos_measurement_chans = {'x126_sn','x127n_sn','x127c_sn','x128n_sn',...
    'x128c_sn','x129n_sn','x129c_sn','x130n_sn','x130c_sn','x131_sn','x131C_sn'};

% Identifiers for the different search classes so that the protein and
% phospho data can be separated
prot_search_class_Col = 'Class';
prot_search_class = 'A';
phos_search_class_Col = 'class';
phos_search_class = 'B';

% Names for two output files 
prot_peps_analyze_file = [dir_str_p1,filesep,date_out,'_',exp_str,'_2018protein_2018ProteinPhosMap.csv'];
phos_peps_analyze_file = [dir_str_p1,filesep,date_out,'_',exp_str,'_2018phos_2018ProteinPhosMap.csv'];

% The number of TMT channels used for the experiment 
data_col_num = 11;

	B. Second … The Python3 matching code. 

	C. Third … “exp_str”_PreOcc_Script_Part2.mlx

		A. DEFINE THE INPUTS: 
			- “exp_str” — the identifier for the experiment — THIS MUST BE WHAT IS WAS FOR PART 1 !!
			- “date_part1” - the date string for Part1 in “yymmdd” format. 
			- “”xen_human_matched_filtered_file” - The name of the file that has the filtered xenopus - human matches cross-reference to the database for low throughput peer-reviewed references 
			
			** The input for the oocyte maturation project is below. 

		B. Press Run. After is finishes check that a new folder was created: “date_str”_”exp_str”_PreOcc_Files_Part2. There should be 1 output file inside of it.   

% Experiment Identification String
exp_str = 'OocyteMat';

% The date string for the date that the desired Part 1 files were generated
% on. 
date_part1 = ‘’; % ("yymmdd")

% The input file that has the human <-> xenopus phospho-residue matches
% filtered by motif score and cross referenced to the human database. 

xen_hum_match_filtered_file = 'OocyteMat_Filtered0p75_LTinfo_01282022.csv';
			
	D. Fourth… “exp_str”_PreOcc_Script_Part3.mlx 

		A. Define the inputs: 

			- “exp_str” — the identifier for the experiment — THIS MUST BE WHAT IS WAS FOR PARTs 1 and 2!!
			- “reference_dir” and “reference_filename” : 
The Xenopus protein model reference and the folder that this reference is in. The reference file and folder name should be the same as they were in the Part1 script. 
			- “date_part1” — The date_str for the Part1 output files. 
			- “date_part2” — The date_str for the Part2 output files. 
			- “data_col_num” — the number of tot channels for the experiment
			
			** The input for the oocyte maturation project is below. 
	
		B. Press Run. After is finishes check that a new folder was created: “date_str”_”exp_str”_PreOcc_Files_Part3. There should be 4 files and 1 mat structure inside.

				-> “date_str”_”exp_str”_MatchedPhosphoSets.csv  is the input file for occupancy estimation. 

				-> “date_str”_”exp_str”_PhosDetrendedOnly_Info_Data.csv is the relative phospho- data that has matching protein data and therefor has been de-trended by the protein data 

				-> “date_str”_”exp_str”_ProteinFractionSignal.csv is the protein data where the peptide values have been summed and the output is the fraction signal. 
		
% Experiment Identification String
exp_str = 'OocyteMat';

% Location and filename for the Xenopus protein model reference 
reference_dir = 'OocyteMat_InputFiles_PreOccupancy';
reference_filename = 'XENLA_9.2_Xenbase_plusMT.pep.fa';

% The date string for the date that the desired Parts 1 and 2 files were generated
% on. 
date_part1 = '220127'; % ("yymmdd")
date_part2 = '220128'; % ("yymmdd")

% The number of TMT channels for the experiment. 
data_col_num = 11;

5) Now we want to estimate occupancy! There are two choices about how to do this.
 
	1) Do locally on your own computer. 
	2) Upload the functions and input files and run on a server. 

This read me contains information about both of these options. The functions and input files are the same for both. However, for option 1) there is a MATLAB live script and for 2) there is a Bash file. Some of the inputs are different depending on which option is chosen. Also, the bash file is written with inputs for a SLURM scheduler. If you are not using a SLURM scheduler then some of these might need to be changed.  

	OPTION 1) “exp_str”_Occ_Script.mlx 

		DEFINE INPUTS
			-  “exp_str” — the identifier for the experiment — THIS MUST BE WHAT IS WAS FOR PARTs 1, 2, and 3!!
			- “date_part3” — The date_str for the Part3 output files. 
			- “output_folder” — The suffix for the output folder. This folder name needs to be different than any current folder that is in the file path, otherwise there will be an error. 
			- “x_axis_label” - The x_axis_label for the experiment
			- “x_axis_points” - The x_axis tick label for the experiment as a comma separated string. 
			- “confidence_interval” - This value * 100 is is the confidence value for the boostrapped confidence intervals
			- “iterations” - The number of bootstrap iterations. Ideally this number is 10,000. However, if this number is chosen and run locally it will take a long time. Often I use 100 if I want to test how a change to the code is working. 
			- “norm_2_protein_FLAG” - Flag for whether you want to normalize to protein. This value should be 1 because that means you are normalizing to protein. If you have very good reason to expect that protein values are not changing you could set the flag to 0. 
			- “fit_plot_FLAG” - Flag for whether you want to plot the occupancy estimation visualizations. These visualizations allow you to see how occupancy estimation works and gives you an intuitive sense of why the confidence intervals look like they do for different data points. However, they increase the time for the code for each phospho-set, so it is recommended to set this flag to 1 only if you trying to get a sense of how estimating occupancy works. 
			- “occ_plot_FLAG” - Flag for whether you want to plot the plot the results of occupancy estimation for each phospho-set. These plots are very useful to identifying phospho-sets that have interesting trends and how tight confidence intervals. Since, in reality, the number of phospho-sets that meets this criteria is low, it is often reasonable to scan these files and pull out the interesting ones. 
			- “cluster_INFO ” - Parameter for whether you are parallelizing  or for how many cores you are using to estimate occupancy locally. 
				0 = no parallelization
				1 = yes parallelize, use all the cores on your computer
				> 1 < 128 = yes parallelize, use a specified number of cores on computer 

			-> Examples of these parameters for oocyte occupancy estimation are below. 
	
	RUN the script and wait for the output to be generated. 
		
		The outputs will be a folder titled: “output_folder”_Output_”date_str”
		This folder will have 2 subfolders (depending on the plot FLAGS and two output files)
			
			If fit_plot_FLAG = 1, there is a folder titled “Fit_plots”
			If occ_plot_FLAG = 1, there is a folder titled “Occ_plots”
			File 1 - ending with “_Matched_PhosphoSets.csv” which has the occupancy results. 
			File 2 - ending with “_Matched_PhosphoSetsexit_codes” which has the exit codes for the occupancy estimation (was occupancy estimated or not, and how many forms were in each phospho-set etc) 

	OPTION 2) ESTIMATE OCCUPANCY ON THE CLUSTER 

	A. Create a new folder that will be uploaded  
		Into this folder put copied versions of 
			- The “_Matched_Phosphosets.csv” file that is in the output of Part 3 folder 
			- Copies of all of the functions inside the “Occ_Functions” folder. Put in all the functions and not the folder. 
			- A bash script (.sh) that will be used to call the functions to estimation occupancy on the server. And example of this bash script is below. 
	
	B. Upload this folder to O2 and execute it. Then download the output folder which will have the same contents described above for running occupancy locally. 
	
		Here are example commands to use to run occupancy on the o2 server using my o2 username which is esv4:

			- Make sure your local directory is the one that has the folder that you are uploading in it. 
			- Upload the contents: scp -r OocyteMat_OccEstServer_220126_v2/ esv4@transfer.rc.hms.harvard.edu:/n/groups/kirschner/evi

			- Sign into the server: ssh esv4@o2.hms.harvard.edu 
				-> verify credentials

			- Navigate to the folder that was just uploaded: cd /n/groups/kirschner/evi/OocyteMat_OccEstServer_220126_v2/

			- Make the .sh script executable: chmod +x 220126_OocyteMat.sh 

			- Start occupancy: sbatch 220126_OocyteMat.sh 220126_OocyteMat_Matched_PhosphoSets.csv 10000 OocyteMat_Occ

			- Check progress: squeue -u esv4

			- Download the output folder:  scp -r esv4@transfer.rc.hms.harvard.edu:/n/groups/kirschner/evi/OocyteMat_OccEstServer_220126_v2/OocyteMat_Occ_10000_iterations_Output_20220126/ ./

~~~~~~~~~~~~~~~ Example Inputs for “occ_est”_Occ_Script.mlx ~~~~~~

% Experiment Identification String
exp_str = 'OocyteMat';

% The date string for the date that the desired Parts 3 file was generated on. 
date_part3 = '220128'; % ("yymmdd")

% The suffix for the output folder; 
output_folder = 'occ_output';

% The x-axis label for the experiment. 
x_axis_label = 'Hours Post Progesterone';

% The x-axis tick lables for the experiment. 
x_axis_points = '0,1,2,3,4,5,6,7,8,ppt-1,ppt-2';

% confidence_interval * 100 is the confidence value for the boostrapped confidence intervals
confidence_interval = 0.9; 

% number of bootstrap iteractions for occupancy estimation
iterations = 100; 

% Flag for whether you want to normalize to protein. 
%   1 = normalize to protein, 0 = do not normalize to protein
norm_2_protein_FLAG = 1; 

% Flag for whether you want to plot the occupancy estimation
% visualizations. 
%   1 = plot, 0 = do not plot
fit_plot_FLAG = 0; 

% Flag for whether you want to plot the occupancy estimation results. 
%   1 = plot, 0 = do not plot 
occ_plot_FLAG = 1; 

% Parameter for whether you are parallelizating or for how many cores you are estimating occupancy with. 
cluster_INFO = 0; 
    % 0 = No parallelization 
    % 1 = Yes parallelize, running locally, use all the cores on your computer. 
    % > 1 < 128 = Yes parallelize, running locally, using a specified number of cores > 1 < 128
    % > 128 = Yes paralellize, run on a cluster [THIS IS NOT AN OPTION FOR
    % THIS SCRIPT BECAUSE IT IS FOR RUNNING LOCALLY ] 

~~~~~~~~~~~~~~~ Example .sh file for estimating occupancy on the server ~~~~~

#!/bin/sh
#SBATCH -p priority
#SBATCH -c 1
#SBATCH -t 1-0
#SBATCH --mem=20G
#SBATCH -o %j.out
module load matlab/2020a
input=$1
iterations=$2
data=$3
matlab -r "tic;addpath('/n/groups/kirschner/evi/220128_OocyteMat_OccEstServer_v1/');modocc('/n/groups/kirschner/evi/220128_OocyteMat_OccEstServer_v1/"$input"',0.90,$iterations,1,0,1,'"$data"_"$iterations"_iterations',129,'Hours Post Progesterone','0,1,2,3,4,5,6,7,8,ppt-1,ppt-2');toc;quit"
