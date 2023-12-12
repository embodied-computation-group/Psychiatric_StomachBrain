function timeseries_preprocessBOLD(subj_idx,cfgMain)

%{

Takes the output of timeseries_prepare_import2matlab that is stored in 
each subject timeseries\data folder 
Y:\Subjects\Subject13\Timeseries\MRItimeseries\fMRItimeseries_S13_kw3.mat

and preprocess it so it can be further analyzed
Preprocessing steps, implemented in fieldtrip, includes remove polinomial trends of second degree,
bandpass filter between 0.01 and 0.1 Hz and standarize to Z units
the output is stored in the tiemseries folder of each subject
Y:\Subjects\Subject13\Timeseries\MRItimeseries\BOLDFULLBAND_S13_kw3


% IR commented on the 12/09/2016
% recommented 28/06/2017

%}

rootDir = global_path2root_folder

filename_bold_input = [global_path2root,'subject',sprintf('%.4d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries',filesep,cfgMain.task,filesep,'BOLDMatlab_rest_fmriprep_smooth_',num2str(cfgMain.kernelWidth),'.mat']

load(filename_bold_input)

filename_bold_output =[global_path2root,'subject',sprintf('%.4d',subj_idx),filesep,'Timeseries',filesep,...
    'MRItimeseries',filesep,cfgMain.task,filesep,'BOLD_filtered_fullband'];

plotDir = strcat(rootDir,'subjects',filesep,'Subject',sprintf('%.4d',subj_idx),filesep,'Plots',filesep);
plotFilename = strcat(plotDir,'S_',sprintf('%.4d',subj_idx),'_',cfgMain.task,'_preprocessBOLD');


BOLDtimeseries.BOLDtsNOPOLI = ft_preproc_polyremoval (BOLDtimeseries.trialVector, 2);
BOLDtimeseries.BOLD_filtered = ft_preproc_bandpassfilter (BOLDtimeseries.BOLDtsNOPOLI,cfgMain.TR,[0.01 0.1]);
BOLD_filtered_zscored = ft_preproc_standardize (BOLDtimeseries.BOLD_filtered);


save(filename_bold_output,'BOLD_filtered_zscored','-v7.3')

%% Sanity check

if cfgMain.savePlots == 1
    
    if cfgMain.plotFigures == 0;
        SanityPlot = figure('visible','off');
    else
        SanityPlot = figure('visible','on');
    end
    

    % TO CHANGE ONCE WE CHANGE THE VOXEL RESOLUTION TO 3MM
%     voxelCoordinates = sub2ind([109,128,108],64,64,82);% for 1mm vmp
    voxelCoordinates = sub2ind([53,63,46],11,30,37);
%[53,63,46],11,30,37)
% voxelCoordinates = sub2ind([79,95,79],9,45,53); % voxel in somatomotor cortex
    
    subplot(4,1,1) % 
    plot(BOLDtimeseries.trialVector(voxelCoordinates,:),'LineWidth',4)
    title(['S',sprintf('%.2d',subj_idx),32,'raw'],'fontsize',18)
    subplot(4,1,2)
    plot(BOLDtimeseries.BOLDtsNOPOLI(voxelCoordinates,:),'LineWidth',4)
    title(['S',sprintf('%.2d',subj_idx),32,'Polinomial removal'],'fontsize',18)
    subplot(4,1,3)
    plot(BOLDtimeseries.BOLD_filtered(voxelCoordinates,:),'LineWidth',4)
    title(['S',sprintf('%.2d',subj_idx),32,'Bandpass filtered between 0.01 and 0.1'],'fontsize',18)
    subplot(4,1,4)
    plot(BOLD_filtered_zscored(voxelCoordinates,:),'LineWidth',4)
    title('zscored','fontsize',18)
    
    
    set(gcf,'units','normalized','outerposition',[0 0 1 1])
    set(gcf, 'PaperPositionMode', 'auto');
    
    print ('-dpng', '-painters', eval('plotFilename'))
    print ('-depsc2', '-painters', eval('plotFilename'))
    saveas(SanityPlot,strcat(plotFilename,'.fig'))
end

end