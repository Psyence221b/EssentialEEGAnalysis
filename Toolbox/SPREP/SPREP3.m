% Preprocessing Semi-automatic One shot tool
% made by SIU
%% load data .cnt file
clc;
% In here, you can load .cnt file in eeglab GUI.
% eeglab;
% setting option not to use single precision.('option_single')
% pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1,...
%     'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
%     'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
%     'option_checkversion', 1, 'option_chat', 0);

%% parameter
HP_filterHz = 0.2; %(0.1 cut-off)

%% Get name and subject number
% divide current EEG file name for search physiology file which has same
% name
[pathstr, name, ext] = fileparts(EEG.comments);
% divide path(folder)
pathsplit = strsplit(pathstr, '\');
physioname = [name,'.mat'];
% current folder location
Rootpath = pwd;

% Start log
tic;
cd(pathsplit{1,4}); mkdir('Log'); cd('Log');
logFile = fopen([pathsplit{1,4} '_Preprocess_log.txt'], 'at+');
cd(Rootpath);
fprintf(logFile, ['\n===================================== \n']);
fprintf(logFile, [' START: ' datestr(now) '\n']);
fprintf(logFile, [' FILE : ' name ' \n']);
fprintf(logFile, [' ¡Ø This file is filter at 0.1 high-pass cutoff frequency \n \n']);

%% load channel location file
%load 64ch + 4 monopolar channels data and add ECG channels
%EEG=pop_chanedit(EEG, 'lookup','C:\Users\pch\Documents\MATLAB\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc','changefield',{60 'labels' 'OI1'},'changefield',{64 'labels' 'OI2'},'lookup','C:\Users\pch\Documents\MATLAB\eeglab14_1_1b\\plugins\\dipfit2.3\\standard_BEM\\elec\\standard_1005.elc','append',68,'changefield',{69 'labels' 'ECG'});
EEG=pop_chanedit(EEG, 'lookup',[pwd '\64ch_4chmonopolar.ced'],'append',68,'changefield',{69 'labels' 'ECG'});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG.nbchan = 69;
EEG.chanlocs(69).labels = 'ECG'; EEG.chanlocs(69).urchan = 69;
EEG.chaninfo.nodatchans = [];
fprintf(logFile, ['%.2f - Channel location set \n'], toc);

%% Merge physiology data file
% change path tempollary
cd('..'); cd('physiology'); cd(pathsplit{1,4});
% load physiology file(.mat)
load(physioname);
% back to previous folder path
cd(Rootpath); 

% calculate MP150 physiology data
data(:,11) = data(:,3)/5*2^7 + data(:,4)/5*2^6 + data(:,5)/5*2^5 + data(:,6)/5*2^4 + data(:,7)/5*2^3 + data(:,8)/5*2^2 + data(:,9)/5*2^1 + data(:,10)/5*2^0;
line_trigger = data(:,11)';
leng = length(line_trigger);

% find first trigger time
[~, ECG_ontrigger] = find(line_trigger~=0);
ECG_firsttrigger = ECG_ontrigger(1,1);

% find first EEG data trigger time
EEG_firsttrigger = EEG.event(1).latency;
leng_EEG = size(EEG.data,2) - EEG_firsttrigger;

% Add physiology data to EEG structure
% add MP150 data to EEG structure
if (ECG_firsttrigger+leng_EEG) <= size(data,1)
    %data_RSP = [zeros(1,Loc_EEGtrigstart-1) data(Loc_firsttrigger:(Loc_firsttrigger+leng_EEG),1)'];
    data_ECG = [zeros(1,EEG_firsttrigger-1) data(ECG_firsttrigger:(ECG_firsttrigger+leng_EEG),2)'];
    fprintf(logFile, ['%.2f - Merge ECG data \n'], toc);
else
    % codes for happening when physiology data is shorter than EEG data
    data_ECG = [zeros(1,EEG_firsttrigger-1) data(ECG_firsttrigger:end,2)'];
    data_ECG = [data_ECG zeros(1,size(EEG.data,2)-size(data_ECG,2))];
    fprintf(logFile, ['%.2f - Merge ECG data (¡Ø ECG is short)\n'], toc);
end

EEG.data(69,:) = data_ECG*200;
EEG = eeg_checkset( EEG );

clear data isi isi_units labels;

%% check trigger number and adjust odd trigger(trigger error)
%load CIT excel data file
cd('..'); cd('CIT'); cd('Excel'); cd(pathsplit{1,4});
if name(end) == '1'
    CIT_name = ['2017_CIT_left_no-',pathsplit{1,4},'-1.xlsx'];
elseif name(end) == '2'
    CIT_name = ['2017_CIT_right_no-',pathsplit{1,4},'-1.xlsx'];
end
CIT_data = xlsread(CIT_name);
cd('..'); cd('..'); cd('..'); cd('EEG');

% adjust EEG.event in another variable temporally
EEGtrigger_temp = [EEG.event(1:end).type];
EEGtrig_loc_start = find(EEGtrigger_temp == 1);
EEGtrig_loc_end = EEGtrig_loc_start - 1;
EEGtrig_loc_end(1) = []; EEGtrig_loc_end = [EEGtrig_loc_end size(EEGtrigger_temp,2)];
EEGtrig_num = EEGtrig_loc_end - EEGtrig_loc_start + 1;
% reshaping trigger vector to matrix form in the 'EEGtrigger '
for i = 1:size(EEGtrig_loc_start,2)
    for j = 1:EEGtrig_num(i)
    EEGtrigger{i,j} = EEGtrigger_temp(EEGtrig_loc_start(i)+(j-1));
    end
end
% change trigger code with behavior data file
for i = 1:size(EEGtrigger,1)
    EEGtrigger{i,2} = CIT_data(i,24); EEGtrigger{i,3} = 3;
end
% change EEGLAB structure with the variable 'EEGtrigger'
for i = 1:size(EEGtrig_loc_start,2)
    for j = 1:EEGtrig_num(i)
    EEG.event(EEGtrig_loc_start(i)+(j-1)).type = EEGtrigger{i,j};
    end
end
EEG = eeg_checkset( EEG );

clear pathstr ext leng leng_EEG line_trigger;
clear Loc_EEGtrigstart Loc_firsttrigger Loc_ontrigger;
clear CIT_data CIT_name EEGtrig_loc_end EEGtrig_loc_start EEGtrig_num EEGtrigger_temp;
fprintf(logFile, ['%.2f - Correct Triggers(Events)  \n'], toc);

%% remove M1, M2 channles and recheck EEGLAB
%remove M1, M2 channels
EEG = pop_select( EEG,'nochannel',{'M1' 'M2'});
%recheck EEGLAB
EEG = eeg_checkset( EEG );
fprintf(logFile, ['%.2f - Remove M1, M2 channels  \n'], toc);

%% Change sampling rate to 250Hz
% It's for speed up to calculation
EEG = pop_resample( EEG, 250);
EEG = eeg_checkset( EEG );
fprintf(logFile, ['%.2f - Downsampling to 250Hz  \n'], toc);

%% High-pass filter
EEG = pop_eegfiltnew(EEG, HP_filterHz, 0);
EEG.setname = [name '__0.1hz'];
fprintf(logFile, ['%.2f - Hi-passed 0.2Hz(0.1 cut-off)  \n'], toc);

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

%% end message
eeglab redraw;
fprintf(logFile, ['\n First Part Preprocessing is done \n \n']);
fclose(logFile); %stop logging

disp('======================================== ');
disp('    First Part Preprocessing is done');
disp('======================================== ');
%%
% save each file, and merge them together to one CIT file