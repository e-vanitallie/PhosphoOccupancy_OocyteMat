function plotOcc_wAll_Info(DataVars, print_path, occ_x_label, occ_x_ticks, ExportTable)
%{
This function substitutes the occupancy plotting protocol at the end of modocc.m (see previous versions of that file).
There is some duplication of code, because as of initial commit, extraction and isolation of the SVD plotting code is borrowing
some variable initializations.
%}

% combine all sets of labels
% Data_labels_combined = [labels_opt,labels_hi,labels_lo, DataVars];

% put the occupancy, hi, low, and ModsetData (parent + phospho, NO
% NORM protein line) into CombinedDataTable
%CombinedDataTable = array2table(round(Combined_data, 2),'VariableNames',Data_labels_combined);
%prot_table = array2table(round(prot_stored, 2),'VariableNames',labels_prot);
%labeltable = ImportData(NoProteinIndex,LabelVars);

% also, want to return plots for ONLY when a hyperplane is fit
% We decided that not_hyperplane_flag is not meaningful as written
% We are going to plot if the total number that have not been corrected
% is at least one.

%    if total_num_uncorrected >= 1

%disp('plotting occupancy output')


labels_opt = cell(1,length(DataVars));
for RenameCounter = 1:length(DataVars)
    labels_opt(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_occ']);
end

labels_hi = cell(1,length(DataVars));
for RenameCounter = 1:length(DataVars)
    labels_hi(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_hi']);
end

labels_lo = cell(1,length(DataVars));
for RenameCounter = 1:length(DataVars)
    labels_lo(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_lo']);
end

labels_prot = cell(1,length(DataVars));
for RenameCounter = 1:length(DataVars)
    labels_prot(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_prot']);
end

labels_orig = cell(1,length(DataVars));
for RenameCounter = 1:length(DataVars)
    labels_orig(1,RenameCounter) = cellstr([DataVars{RenameCounter},'_orig']);
end

% Enter appropriate x axis for plots
% x_axis_data contains 1 through 11
x_axis_data = 1:1:size(DataVars,2);

% Enter colors
colors = lines;

%figure this out later
%plot_fitting_scatter = 1;

%x_axis_labeling = 'Time Intervals (hours)';
x_axis_labeling = occ_x_label;

%ConditionLabels = {'0';'2';'3';'4';'5';'6';'7';'8';'9';'(+pptase) 10';'(+pptase) 11'};
ConditionLabels =  occ_x_ticks;
% ConditionLabels = {};

% table containing redundant list of set ids
set_ids_table = ExportTable.set_id;
site_ids_table = ExportTable.P_Form;

% unique list of set IDs - should be only one here, with rewritten modocc

%SetsOfSitestable = unique(ExportTable.set_id); % since only once input this is just 1

% redundant list of protein IDs
protid_plot = ExportTable.Protein_Reference;
% redundant list of gene symbols
gs_plot = ExportTable.GeneSymbol;
gs_x_plot = ExportTable.GeneSymbol_Xenbase;
gs_m_plot = ExportTable.Phosphosite_Motif;
gs_r_plot = ExportTable.Phosphosite_Position;

hm_gs_plot = ExportTable.PP_Reference;
hm_res_plot = ExportTable.Human_Residue;
hm_m_plot = ExportTable.Human_Motif;
hm_info_plot = ExportTable.MoreInfoFlag;

% is there protein data
prot_plot_FLAG = ExportTable.Match_Protein;

% extract data for occupancy, CIs, ModsetData (parent + phospho), and
% norm protein to separate tables
opt_data_plot = ExportTable{:,labels_opt};
hi_data_plot = ExportTable{:,labels_hi};
lo_data_plot = ExportTable{:,labels_lo};
rel_data_plot = ExportTable{:,DataVars};
orig_data_plot = ExportTable{:,labels_orig};
prot_data_plot = ExportTable{:,labels_prot};

%for plot_counter = 1:size(SetsOfSitestable,1)

% first plot for protein
% then plots for each occupancy
% last the relative trend plots

Set_table_indicies = 1:1:size(set_ids_table,1); %(plot_counter)==set_ids_table;

rel_trends = rel_data_plot(Set_table_indicies,:);
opt_trend = opt_data_plot(Set_table_indicies,:);
hi_trend = hi_data_plot(Set_table_indicies,:);
lo_trend = lo_data_plot(Set_table_indicies,:);
orig_trend = orig_data_plot(Set_table_indicies,:);
prot_trend = prot_data_plot(Set_table_indicies,:);

rel_trends(rel_trends==0)=1E-9;


protid_label = protid_plot(Set_table_indicies);
gs_label = gs_plot(Set_table_indicies);
gs_x_label = gs_x_plot(Set_table_indicies);
gs_m_label = gs_m_plot(Set_table_indicies);
gs_r_label = gs_r_plot(Set_table_indicies);

hm_gs_label = hm_gs_plot(Set_table_indicies); 
hm_res_label = hm_res_plot(Set_table_indicies);
hm_m_label = hm_m_plot(Set_table_indicies);
hm_info_label = hm_info_plot(Set_table_indicies); 

NumForms2 = size(opt_trend,1);

if (NumForms2 > 1 && NumForms2 <= 10) %% this should already have been caught by bootci && (NumForms2 < NumConditions)
    
    %if max((max(hi_trend-lo_trend)))>0
    
    occfig = figure('visible','on');
    
    %sgtitle(['Set ',num2str(SetsOfSitestable(plot_counter)),' ',protid_label{1},', ', gs_label{1}], 'fontsize', 36);
    %sgtitle([protid_label{1},', ', gs_label{1}], 'fontsize', 36);
    
    if (isempty(hm_gs_label{1}))
        title_str = [gs_x_label{1},', ',gs_r_label{NumForms2},', ',gs_m_label{NumForms2},', ',gs_label{1}];
        sgtitle(title_str, 'fontsize', 25);
    else
        human_split = split(hm_gs_label{1},"|");
        
        moreinfo_string = 'No'; 
        for c = 2:NumForms2
            if regexp(hm_info_label{c},'NaN')
            else
                moreinfo_string = 'Yes';
                break
            end
        end
                
                
        title_str = {['Xen: ',gs_x_label{1},', ',gs_r_label{NumForms2},', ',gs_m_label{NumForms2},', ',gs_label{1}],...
            ['Human: ',human_split{2},', ',hm_res_label{2},', ',hm_m_label{2},', More Info: ',moreinfo_string]};
        sgtitle(title_str, 'fontsize', 25);
    end
    
    numrows = ceil((NumForms2 + 2)/ 4);
    
    %% FIRST PLOT THE PROTEIN
    
    subplot(numrows,4,1)
    
    if ( prot_plot_FLAG )
        
        plot_prot_trend = (prot_trend(1,:));
        
        plot(x_axis_data,plot_prot_trend,'color','k','linestyle','-','linewidth',1.5,'marker','.','markersize',10)
        xlabel(x_axis_labeling)
        ylabel('Fraction of total protein signal')
        
        ylim([0,(max(max(plot_prot_trend)))*1.1])
        
    else
    end
    
    set(gca,'fontsize',12)
    
    if isempty(ConditionLabels)==0
        set(gca,'xtick',x_axis_data,'xticklabel',ConditionLabels)
        set(gca,'XTickLabelRotation',45)
    end
    
    title('Protein trend','fontsize',12)
    set(gca, 'FontName', 'Arial')
    
    box off
    set(gca,'Layer', 'Top')
    ax = gca;
    c = ax.PlotBoxAspectRatio;
    ax.PlotBoxAspectRatio = [0.7043    1.0000    0.7043];
    set(gcf,'color','w')
    set(gca, 'FontName', 'Arial')
    %%
    
    legend_names={};
    
    for form_counter = 1:NumForms2
        
        subplot(numrows,4,1+form_counter) % the first plot was protein,and now we count up from there
        
        upper = hi_trend((form_counter),:);
        lower = lo_trend((form_counter),:);
        optimal_trend = opt_trend((form_counter),:);
        
        optimal_trend(isnan(optimal_trend))=0;
        
        
        fill([x_axis_data, fliplr(x_axis_data)], [upper, fliplr(lower)], colors(form_counter,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on
        plot(x_axis_data, optimal_trend,'color',colors(form_counter,:),'linestyle','-', 'linewidth',1.5,'Marker','.','markersize',12)
        hold on
        
        %         sitepositionnum = single_Phos_and_NON_phos{site_counter,'SitePosition'};
        %         pretitle2 = single_Phos_and_NON_phos{site_counter,'GeneSymbol'};
        %         text(x_axis_data(end)+0.25,optimal_trend(end)+2,[pretitle2{1,1},', Site ',num2str(sitepositionnum)],'fontsize',12,'color','k')
        xlabel(x_axis_labeling)
        ylabel('% Occupancy   ')
        ylim([0,100])
        set(gca,'fontsize',12)
        
        if (form_counter==1)
            legend_names{1}='Un-phosphorylated';
            title('Un-phosphorylated   ','fontsize',12)
        else
            
            title(['Residue(s) ',gs_r_label{form_counter},'   '],'fontsize',12)   % this should use data in phospho-site
            
            legend_names{form_counter} = ['Residue(s): ',gs_r_label{form_counter}];
            
        end
        
        set(gca, 'FontName', 'Arial')
        
        if isempty(ConditionLabels)==0
            set(gca,'xtick',x_axis_data,'xticklabel',ConditionLabels)
            set(gca,'XTickLabelRotation',45)
        end
        
        box off
        
        set(gca,'Layer', 'Top')
        %x_axis_data = 0:2:18;
        ax = gca;
        c = ax.PlotBoxAspectRatio;
        ax.PlotBoxAspectRatio = [0.7043    1.0000    0.7043];
        set(gcf,'color','w')
        set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
        
    end % end form_counter loop
    
    
     %%% This is code to plot the relative trends. This should
    %%% be the last THING plotted because the least informative
    
    subplot(numrows,4,NumForms2+2)
    hold on 
    plot_rel_trends = (rel_trends);
    colors_rel = [0.7 0.7 0.7; 0.5 0.5 0.5; 0.3 0.3 0.3];
    colors_rel = [colors_rel; zeros(7,3)];
    
    for form_counter = 1:NumForms2
        
        %plotModform = plot_rel_trends(form_counter,:);
        plotModform = plot_rel_trends(form_counter,:)./sum(plot_rel_trends(form_counter,:));
        plotOrigform = orig_trend(form_counter,:);
        
        plot(x_axis_data, plotModform,'color',colors_rel(form_counter,:),'linestyle','-','linewidth',1.5,'marker','.','markersize',10,...
            'DisplayName',[legend_names{form_counter},' - Detrended'])
         plot(x_axis_data, plotOrigform,'color',colors_rel(form_counter,:),'linestyle','--','linewidth',1,'marker','.','markersize',3,...
            'DisplayName',legend_names{form_counter})
       
        
        if isempty(ConditionLabels)==0
            % gca - current axes or chart for current figure
            set(gca,'xtick',x_axis_data,'xticklabel',ConditionLabels)
            set(gca,'XTickLabelRotation',45)
        end
        
        %         sitepositionnum = single_Phos_and_NON_phos{site_counter,'SitePosition'};
        %         pretitle2 = single_Phos_and_NON_phos{site_counter,'GeneSymbol'};
        
        % short name for black is k
        %text(x_axis_data(end)+0.2,plotModform(end),num2str(form_counter-1),'fontsize',12,'color','k')
    end
    
    xlabel(x_axis_labeling)
    ylabel('Fraction Total Signal')
    
%     if min(min(plot_rel_trends)) < -1 || max(max(plot_rel_trends)) > 1
%         ylim([(min(min(plot_rel_trends)))-0.5,(max(max(plot_rel_trends)))+0.5])
%     else
%         ylim([-1,1])
%     end
    
    set(gca,'fontsize',12)
   % legend('Location','northwestoutside');
     legend('Position', [0.1 0.8 0.05 0.1]);
   % legend('Position', [0.7 0.8 0.7 0.1]);
    legend show
    title('Relative Trends','fontsize',12)
    set(gca, 'FontName', 'Arial')
    
    % box off removes the box outline around the current axes by setting their Box property to 'off'
    % https://www.mathworks.com/help/matlab/ref/box.html?searchHighlight=box%20off&s_tid=doc_srchtitle
    box off
    
    set(gca,'Layer', 'Top')
    ax = gca;
    c = ax.PlotBoxAspectRatio;
    ax.PlotBoxAspectRatio = [0.7043    1.0000    0.7043];
    % gcf - current figure handle
    set(gcf,'color','w')
    
    %%%%%%%%% finishing code for the relative trend %%%%%%%%
    
    setname = set_ids_table(1);
    prot_ref_split = split(protid_label{1},"|");
    
    % comment this line if you do not want the PNG per usual modocc operation
    print(occfig,'-dpng',[print_path,'Occupancy_plots',filesep,'Set_',num2str(setname),'_',prot_ref_split{3},'_',gs_label{1}, '_', num2str(NumForms2),'_forms_occ.png'])
   % close(occfig)
    % uncomment this line to save as fig - which you can open up in Matlab and edit if you wish
    %savefig(occfig, [print_path,'Set_',num2str(setname),'_',protid_label{1},'_',gs_label{1},'_occ.fig'])
    %print(occfig,'-depsc',['Set_',num2str(setname),'_',protid_label{1},'_occ_',exp_name])
    
    %   end % end NumForms2 if
    
   
    
    %    end % end for the if statement for the not_hyperplane_flag check
    % function end
end
end

