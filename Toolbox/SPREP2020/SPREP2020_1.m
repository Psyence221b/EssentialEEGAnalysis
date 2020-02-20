% Semi-automatic EEG preprocessing tool written by SIU
% Revised at Feb, 2020
% -------------------------------------------
% - change reference from linked-earlobes to average earlobes
% - using EEGLAB basic fucntion ICLable

addpath(pwd);

%% Parameter Settings
% Sampling rate
p.SampleRate = 250;

% EEG filtering
p.LP_EEGfilter = 89; %(100 cut-off)
p.HP_EEGfilter = 0.2; %(0.1 cut-off)

% EOG filtering
p.LP_EOGfilter = 89; %(100 cut-off)
p.HP_EOGfilter = 1; %(0.5 cut-off)

% EKG filtering
p.LP_EKGfilter = 89; %(100 cut-off)
p.HP_EKGfilter = 1; %(0.5 cut-off)

%% Load data .cnt file
% In here, please import your raw EEG file
clc;

%% Get name and subject number
% name
[pathstr, name, ext] = fileparts(EEG.comments);
% divide path(folder)
pathsplit = strsplit(pathstr, '\');
% current folder location
Rootpath = pwd;

clear pathstr

%% Start log
tic;
cd(pathsplit{1,end}); 
foldername = 'Log';
mkdir(foldername); cd(foldername);
logFile = fopen([pathsplit{1,end} '_Preprocess_log.txt'], 'at+');
cd(Rootpath);
fprintf(logFile, ['\n===================================== \n']);
fprintf(logFile, [' START: ' datestr(now) '\n']);
fprintf(logFile, [' FILE : ' name ' \n \n']);

%% load channel location file
% load 62ch location file (.ced)
EEG=pop_chanedit(EEG, 'lookup',[pwd '\Neuroscan_62ch.ced']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

fprintf(logFile, ['%.2f - Channel location set \n'], toc);

%% Change sampling rate to 250Hz
% It's for speed up to calculation
EEG = pop_resample( EEG, p.SampleRate);
EEG = eeg_checkset( EEG );
fprintf(logFile, ['%.2f - Downsampling to %dHz  \n'], toc, p.SampleRate);

%% Re-reference to average earlobes except 3 bipolar channels
% Channel numbers 33, 43 are earlobe1, earlobe2 each
% Channels 65~67 are HEOG, VEOG, EKG (bipolar channels)
EEG = pop_reref( EEG, [33 43] ,'exclude',[65:67] );

%recheck EEGLAB
EEG = eeg_checkset( EEG );
fprintf(logFile, ['%.2f - Re-reference to average earlobes  \n'], toc);

%% EEG band-pass filtering
EEG = pop_eegfiltnew(EEG, 'hicutoff',p.LP_EEGfilter,'channels',{'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F7' 'F5' 'F3' 'F1' 'FZ' 'F2' 'F4' 'F6' 'F8' 'FT7' 'FC5' 'FC3' 'FC1' 'FCZ' 'FC2' 'FC4' 'FC6' 'FT8' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'TP8' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO5' 'PO3' 'POZ' 'PO4' 'PO6' 'PO8' 'CB1' 'O1' 'OZ' 'O2' 'CB2'});
EEG = pop_eegfiltnew(EEG, 'locutoff',p.HP_EEGfilter,'channels',{'FP1' 'FPZ' 'FP2' 'AF3' 'AF4' 'F7' 'F5' 'F3' 'F1' 'FZ' 'F2' 'F4' 'F6' 'F8' 'FT7' 'FC5' 'FC3' 'FC1' 'FCZ' 'FC2' 'FC4' 'FC6' 'FT8' 'T7' 'C5' 'C3' 'C1' 'CZ' 'C2' 'C4' 'C6' 'T8' 'TP7' 'CP5' 'CP3' 'CP1' 'CPZ' 'CP2' 'CP4' 'CP6' 'TP8' 'P7' 'P5' 'P3' 'P1' 'PZ' 'P2' 'P4' 'P6' 'P8' 'PO7' 'PO5' 'PO3' 'POZ' 'PO4' 'PO6' 'PO8' 'CB1' 'O1' 'OZ' 'O2' 'CB2'});
EEG.setname = [name '_0.1Hz'];
fprintf(logFile, ['%.2f - EEG Band-passed at %d Hz ~ %d Hz (Cut-off : 0.1 ~ 100)  \n'], toc, p.HP_EEGfilter, p.LP_EEGfilter);

%% EOG band-pass filtering
EEG = pop_eegfiltnew(EEG, 'hicutoff',p.LP_EOGfilter,'channels',{'HEO' 'VEO'});
EEG = pop_eegfiltnew(EEG, 'locutoff',p.HP_EOGfilter,'channels',{'HEO' 'VEO'});
fprintf(logFile, ['%.2f - EOG Band-passed at %d Hz ~ %d Hz (Cut-off : 0.5 ~ 100)  \n'], toc, p.HP_EOGfilter, p.LP_EOGfilter);

%% EKG band-pass filtering
EEG = pop_eegfiltnew(EEG, 'hicutoff',p.LP_EKGfilter,'channels',{'EKG'});
EEG = pop_eegfiltnew(EEG, 'locutoff',p.HP_EKGfilter,'channels',{'EKG'});
fprintf(logFile, ['%.2f - EKG Band-passed at %d Hz ~ %d Hz (Cut-off : 0.5 ~ 100)  \n'], toc, p.HP_EKGfilter, p.LP_EKGfilter);

EEG = eeg_checkset( EEG ); eeglab redraw;

%% Remove 60Hz Line Noise
% cleanLineNoise
% for parallel pool (parpool) function working to execute PREP, 
distcomp.feature( 'LocalUseMpiexec', false );

% line noise correction using 'cleanLinenoise' in PREP
lineNoiseIn = struct('lineNoiseChannels',[1:EEG.nbchan],'lineFrequencies',[60,120,180,240],...
    'fPassBand',[0 EEG.srate/2],'fScanBandWidth',2,'Fs',EEG.srate,'maximumIterations',15,'p',0.01,...
    'pad',0,'taperBandWidth',2,'taperWindowSize',5,'taperWindowStep',1,'tau',100);
[EEG, lineNoiseOut] = cleanLineNoise(EEG, lineNoiseIn);

EEG = eeg_checkset( EEG );
clear lineNoiseIn lineNoiseOut i j;
fprintf(logFile, ['%.2f - Remove 60Hz line noise   \n'], toc);

%% Trigger number change
ChangeTrigger_2020IIF

%% Save 'pre-' file
cd(pathsplit{1,end}); 
savename = [name '_pre.set']
EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
cd(Rootpath);

%% end message
eeglab redraw;
fprintf(logFile, ['\n First Part Preprocessing is done \n \n']);
fclose(logFile);

disp('======================================== ');
disp('    First Part Preprocessing is done');
disp('======================================== ');

clear ext foldername HP_EEGfilter LP_EEGfilter HP_EOGfilter LP_EOGfilter HP_EKGfilter LP_EKGfilter 