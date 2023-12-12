function  timeseries_csfSignal_obtainAndRegress(subj_idx,cfgMain)

%{

This function retrieves the timeseries of a 3x3 cube of voxels located in
the CSF and regress it out the of all the brain voxels, 
the timeseries, the residuals of the regression and the betas coefficient are all stores in the timeseries
folder of each subject.
The rest of PHYSIENS analysis are performed on the residuals of this
regression, stored in the timeseries/data/ folder with the
csfRegressionResiduals_FB_S_SUBJECTNUMER name
The timeseries scripts have the suffix _regression to reflect. 

Output
Bold timeseries csf_regressed
Y:\Subjects\Subject13\Timeseries\MRItimeseries\csfRegressionResiduals_FB_S_13_kw3

IR commented the 12/09/2016
Checked 28/06/2017

%}

%% Input output filenames

% Timeseries
filename_csfr_Residuals_FB = [global_path2root,'subject',sprintf('%.4d',subj_idx),filesep,'Timeseries',filesep,...
    'MRItimeseries',filesep,cfgMain.task,filesep,'BOLD_filtered_fullband_CSFWMREGRESSED']; 

plotDir = strcat(rootDir,'subjects',filesep,'Subject',sprintf('%.4d',subj_idx),filesep,'Plots',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.2d',subj_idx),'_',cfgMain.task,'_CSFregression');

%% New VMP

% load tsv from fmriprep

regressorsfmriprep=readtable('/Volumes/Seagate/Physio_VMP/data/MRI/sub-0019_ses-session1_task-rest_run-1_desc-confounds_timeseries.csv');
%%

insideBrain = tools_getIndexBrain('inside'); % voxels inside brain


filename_bold = [global_path2root,'subject',sprintf('%.4d',subj_idx),filesep,'Timeseries',filesep,...
    'MRItimeseries',filesep,cfgMain.task,filesep,'BOLD_filtered_fullband']; %global_filename(subj_idx,cfgMain,'BOLD_filtered_fullbandFilename'); % Load BOLD timeseries
load(filename_bold)

% REGRESSION

toBeExplained = BOLD_filtered_zscored'; % BOLD timeseries will be the variable to be explained out in the GLM
csf_timeseries = regressorsfmriprep.csf;
wm_timeseries= regressorsfmriprep.white_matter;
regressors=[wm_timeseries,csf_timeseries];


regressorsNOPOLI = ft_preproc_polyremoval (regressors', 2);
regressorsNOPOLI_filtered = ft_preproc_bandpassfilter (regressorsNOPOLI,cfgMain.TR,[0.01 0.1]);
regressors_filtered_zscored = ft_preproc_standardize (regressorsNOPOLI_filtered);

% 
% figure;
% subplot(4,1,1)
% plot(regressors)
% subplot(4,1,2)
% plot(regressorsNOPOLI')
% subplot(4,1,3)
% plot(regressorsNOPOLI_filtered')
% subplot(4,1,4)
% plot(regressors_filtered_zscored')
% 
regressors_filtered_zscored=regressors_filtered_zscored';

% explainingVariables = zscore([csf_timeseries]',[],1); % Variables used to explain the BOLD data
betas_csf_wm = tools_EfficientGLM(toBeExplained,regressors_filtered_zscored); % Obtain the betas indicating how much the predicting variable predicts the data
predictedBOLD = regressors_filtered_zscored(:,1) *betas_csf_wm(1,:)+regressors_filtered_zscored(:,2) *betas_csf_wm(2,:); % What the BOLD timeseries should look like if CSF predicted at 100% accuracy the data
error_csf_wm = toBeExplained - predictedBOLD; % The error is the portion of the data not predicted by the CSF signal
error_csf_wm_z = zscore(error_csf_wm,[],1);
error_csf_wm_z = error_csf_wm_z(:,insideBrain); 
betas_csf_wm = betas_csf_wm(:,insideBrain);


%% Plot regression to check it works

if cfgMain.savePlots == 1

voxelCoordinates = sub2ind([53,63,46],11,30,37);
voxelCoordinates_inside = zeros(53*63*46,1);
voxelCoordinates_inside(voxelCoordinates)=1;
voxelCoordinates_inside = voxelCoordinates_inside(insideBrain);

if cfgMain.plotFigures == 0;
    SanityPlot = figure('visible','off');
else
    SanityPlot = figure('visible','on');
end


hold on
subplot(2,1,1)
plot(regressors,'LineWidth',2)
title(['S',sprintf('%.2d',subj_idx),32,'CSF timeseries'],'fontsize',18)
grid on
subplot(2,1,2)
plot(toBeExplained(:,voxelCoordinates),'r-','LineWidth',4)
hold on
plot(error_csf_wm_z(:,logical(voxelCoordinates_inside)),'b--','LineWidth',3)
grid on
legend ('Before regression','After regression')
title(['S',sprintf('%.2d',subj_idx),32,'EFFects of CSF regressionin right SS'],'fontsize',18)

   
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf, 'PaperPositionMode', 'auto');

print ('-dpng', '-painters', eval('plotFilename'))
print ('-depsc2', '-painters', eval('plotFilename'))
saveas(SanityPlot,strcat(plotFilename,'.fig'))

end

%% SAVE
% save(filename_csfSignal_FB,'csf_timeseries')
save(filename_csfr_Residuals_FB,'error_csf_wm_z')
% save(filename_csfrBetas_FB,'betas_csf_wm')

end