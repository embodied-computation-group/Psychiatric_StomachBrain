function timeseries_prepare_import2matlab(subj_idx,cfgMain)

%{
Construct fMRI timeseries (time,voxels) from swaf images located at each subject folder
and saves them as a structure in the HDD.

Input files
hdr files from each volume
Y:\Subjects\Subject13\fMRI\acquisition1\RestingState\s3wafPHYSIENS_Sujet13-0002-00001-000001-01.hdr

Output 
concatenated MRI timeseries
Y:\Subjects\Subject13\Timeseries\MRItimeseries\fMRItimeseries_S13_kw3.mat

Ignacio Rebollo 16/3/2015
commented 28/06/2017

%}

% kernelWidth = cfgMain.kernelWidth;

% Import

isrest = strcmp(cfgMain.task,'REST')
if isrest
dataDir='/Volumes/Seagate/Physio_VMP/data/MRI/2023/sub-0019/ses-session1/'; %'/Volumes/Seagate/Physio_VMP/data/MRI/'

% dataDir = strcat(global_path2subject(subj_idx,cfgMain.sample),'fMRI',filesep); %navi

% dataDir = strcat(global_path2subject(subj_idx),'fMRI',filesep,'acquisition1',filesep,'RestingState',filesep); %phys

else
% dataDir = strcat(global_path2subject(subj_idx,cfgMain.sample),'taskfMRI',filesep,cfgMain.task);
end

output_filename = [global_path2root,'subject',sprintf('%.4d',subj_idx),filesep,'Timeseries',filesep,'MRItimeseries',filesep,cfgMain.task,filesep,'BOLDMatlab_rest_fmriprep_smooth_',num2str(cfgMain.kernelWidth)]

BOLDtimeseries = []; % "raw" fmri data
BOLDtimeseries.fsample  = 1/cfgMain.TR;  %0.5 hz, ~ TR=2s


% identify the files corresponding to the deisred spatial kernel smoothing
% obtained from the output from mri preprocessing 


%     filename= dir( fullfile( dataDir,'s3_sub-',sprintf('%.4d',subj_idx),'_ses-session1_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii')); %# list all *.hdr files of preprocessed images, this are the 450 volumes
%     filename= dir( fullfile( dataDir,'s3waf*.img')); %# list all *.hdr files of preprocessed images, this are the 450 volumes

filename=[dataDir,'func/', 'sub-',sprintf('%.4d',subj_idx),'_ses-session1_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold_3mm.nii'] %# list all *.hdr files of preprocessed images, this are the 450 volumes]
BOLDtimeseries =ft_read_mri(filename);
dimentions_anatomy=size(BOLDtimeseries.anatomy);
BOLDtimeseries.trialVector=reshape(BOLDtimeseries.anatomy,dimentions_anatomy(1)*dimentions_anatomy(2)*dimentions_anatomy(3),dimentions_anatomy(4));



BOLDtimeseries.time  = [0:cfgMain.TR:cfgMain.nVolumes*cfgMain.TR-1]; % create time axis

save(output_filename,'BOLDtimeseries','-v7.3')


end