% Semi-automatic EEG preprocessing tool written by SIU
% Revised at Feb, 2020
% -------------------------------------------
% - change reference from linked-earlobes to average earlobes
% - using EEGLAB basic fucntion ICLable

addpath(pwd);

%% set parameter
p.TMP_epochrange        = [-0.8 1.5];       % unit(second)
p.Range_microdetect     = [-200 1200];      % ms
p.AMICAthreads          = 8;                % full CPU threads number of your computer
p.AMICAsave             = 1;                % if 1 = save .set file // 0 = no save

p.Range_microdetect     = [-200 1200];      % ms
p.Usemicro              = 0;                % if 1 = Doing microsaccade procedure // 0 = not (If your target is P3, I recommend 0 )
p.Usemicrorange         = 0;                % if 1 = using defined microdetect range // 0 = whole epoch range

p.Numblock              = 6;                % number of blocks
p.Numtrial              = 62;               % number of trials within each block
p.Epochtrigger          = cell(1,p.Numblock*p.Numtrial);
for i = 1:size(p.Epochtrigger,2)
p.Epochtrigger{i}       = [num2str(i) '-sti'];
end

%% Get name and subject number
% divide current EEG file name 
[pathstr, name, ext] = fileparts(EEG.comments);
% divide path(folder)
pathsplit = strsplit(pathstr, '\');
% current folder location
Rootpath = pwd;
clear pathstr ext

%% Keeping original EEG data
% (for interpolation later)
KeepEEG = EEG; 
tic;

%% Start Log
cd(pathsplit{1,end});
%foldername = ['Log_' name];
foldername = 'Log';
mkdir(foldername); cd(foldername);
logFile = fopen([pathsplit{1,end} '_Preprocess_log.txt'], 'at+');
fprintf(logFile, ['\n===================================== \n']);
fprintf(logFile, [' START: ' datestr(now) '\n']);
fprintf(logFile, [' FILE : ' EEG.setname ' \n \n']);

%% Remove Extreme bad channels
% if you have systemical "bad" channels, first remove it.
% ex) really noisy FP1 channel. 
eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Select Extreme bad channels', ...
                 'eloc_file', EEG.chanlocs,'events',EEG.event,'winlength',70,'spacing',150);

EEG = pop_select(EEG);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
EEG.setname = [name '_remove Extreme channel'];
ALLEEG(CURRENTSET).setname = [name '_REJchan1'];


% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('Extreme chan remove', 'Confirm');
eeglab redraw;

% Save 'Remove Extreme bad channels(if exist)'
RMchan1 = []; CHAN_whole = {}; CHAN_remove = {}; REJchan_1 = [];
if KeepEEG.nbchan ~= EEG.nbchan
    % find removed channel information 
    for j = 1:size(KeepEEG.chanlocs,2)
    CHAN_whole{1,j} = KeepEEG.chanlocs(j).labels; end
    for k = 1:size(EEG.chanlocs,2)
    CHAN_remove{1,k} = EEG.chanlocs(k).labels; end
    %find channel name
    RMchan1 = setdiff(CHAN_whole, CHAN_remove,'stable'); 
    for j = 1:size(RMchan1,2)
    REJchan_1(j) = find(strcmp({KeepEEG.chanlocs.labels}, RMchan1{j})==1); end
    % there is a function, eeg_decodechan. remain the code for reference.
    % [~, RMChanlist] = eeg_decodechan(EEGdummyIn.chanlocs, indelec1);
    fprintf(logFile, ['%.2f - Channel removed by Visual inspection: '], toc);
    for ChanList = 1:size(REJchan_1,2)
    fprintf(logFile, '%s(%d)\t', RMchan1{ChanList},REJchan_1(ChanList)); %  \t make space
    end
    fprintf(logFile, '\n');
else
    fprintf(logFile, ['%.2f - Channel removed by Visual inspection: None \n'], toc);
    REJchan_1=[];
end
    % save channel list to mat file
    save('REJchan_1', 'REJchan_1');
    
clear ChanList RMChanlist CHAN_remove CHAN_whole j k

%% Visual Inspection
% In here, plotting EEG data, drag areas where the subject had rested.
% Plot - Channel data(scroll) : Time range to display to 80
pop_eegplot( EEG, 1, 1, 1,[],'winlength',70,'spacing',60);

% make 'Continue button' to wait until VIsual Inspection process is done
SPREP2020_Nextbutton('Visual Inspection', 'Confirm');

% set name current file
EEG.setname = [name '_Visual_Inspection'];
ALLEEG(CURRENTSET).setname = [name '_Visual_Inspection'];
eeglab redraw;

% Save rejected continuous data by visual inspection
REJcont_1 = [];
if ~isempty(TMPREJ)
REJcont_1 = TMPREJ(:,1:2); end
save ('REJcont_1','REJcont_1');
fprintf(logFile, ['%.2f - Continuous data rejection by Visual inspection: REJcont_1.mat \n'], toc);
TMPREJ = [];

%% Find Bad Channels
% find bad channels by using "robust reference" in PREP
EEGdummyIn = EEG; 
referenceIn = struct('referenceChannels',[1:size(EEG.data,1)-3],'evaluationChannels',[1:size(EEG.data,1)-3],...
    'rereference',[1:size(EEG.data,1)-3],'referenceType','robust','channelLocations',EEG.chanlocs,...
    'channelInfo',EEG.chaninfo,'srate',EEG.srate,...
    'robustDeviationThreshold',7,'highFrequencyNoiseThreshold',4.5,'correlationThreshold',0.3,'ransacOff',true);
[EEGdummyOut, referenceOut] = performReference(EEGdummyIn,referenceIn);
disp('   PREP bad channels search is done!');

% ----------------------------------------------
% additional find High-frequency noisy channel
% find channel number F7, F5, or F3
TMP_start = find(strcmp({EEG.chanlocs.labels},'F7')==1);
if isempty(TMP_start), TMP_start = find(strcmp({EEG.chanlocs.labels},'F5')==1);
if isempty(TMP_start), TMP_start = find(strcmp({EEG.chanlocs.labels},'F3')==1); 
if isempty(TMP_start), TMP_start = find(strcmp({EEG.chanlocs.labels},'F1')==1); 
end, end, end

% find channel number P8, P6, or P4
TMP_end = find(strcmp({EEG.chanlocs.labels},'TP8')==1);
if isempty(TMP_end), TMP_end = find(strcmp({EEG.chanlocs.labels},'CP6')==1);
if isempty(TMP_end), TMP_end = find(strcmp({EEG.chanlocs.labels},'CP4')==1);    
if isempty(TMP_end), TMP_end = find(strcmp({EEG.chanlocs.labels},'CP2')==1);    
end, end, end

% find bad channels additionaly using 'Automatic channel rejection'
TMP_list=[];
[TMP_EEG, indelec1]= pop_rejchan(EEG, 'elec',[1:EEG.nbchan-3] ,'threshold',3,'norm','on','measure','spec','freqrange',[0 1] );
[TMP_EEG, indelec2]= pop_rejchan(EEG, 'elec',[TMP_start:TMP_end] ,'threshold',2.5,'norm','on','measure','spec','freqrange',[50 100] );
indelec2 = indelec2 + TMP_start - 1; %adjust channel number because it start in TMP_start(F7)
TMP_list = unique([TMP_list referenceOut.badChannels.all indelec1 indelec2]);
% find removed channels name
[~, RMchan2] = eeg_decodechan(EEG.chanlocs, TMP_list);
clear TMP_EEG TMP_start TMP_end 
disp(' =============================== ');
disp('Find Bad Channels process is done! ');
disp(' =============================== ');

% plotting EEG signal panel to see 'reject' channels
if ~isempty(TMP_list)
    chanlist = [1:EEG.nbchan-3];
    bad_colors = cell(1,length(chanlist)); bad_colors(:) = { 'k' };
    bad_colors(TMP_list) = { 'r' }; bad_colors = bad_colors(end:-1:1);

    tmpcom = [ 'EEGTMP = pop_select(EEG, ''nochannel'', [' num2str(chanlist(TMP_list)) ']);' ];
    tmpcom = [ tmpcom ...
                '[ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
                '   if ~isempty(tmpcom),' ... 
                '     EEG = eegh(LASTCOM, EEG);' ...
                '     eegh(tmpcom);' ...
                '     eeglab(''redraw'');' ...
                '  end; clear EEGTMP tmpcom;' ];

    eegplot(EEG.data(chanlist,:,:), 'srate', EEG.srate, 'title', 'Automatically selected channels to reject', ...
                 'limits', [EEG.xmin EEG.xmax]*1000, 'color', bad_colors(end:-1:1), 'eloc_file', EEG.chanlocs(chanlist), 'command', tmpcom,...
                 'winlength',30,'spacing',40);
  
    % make 'Continue button' to wait until channel rejection process is done
    SPREP2020_Nextbutton('Find Bad channels', 'Confirm');
    
    % set name current file
    EEG.setname = [name '_REJchan2'];
    ALLEEG(CURRENTSET).setname = [name '_REJchan2'];
    eeglab redraw;
    
    % find removed channels number in original data
    for j = 1:size(RMchan2,2)
    REJchan_2(j) = find(strcmp({KeepEEG.chanlocs.labels}, RMchan2{j})==1); end
    % log file automatically removed channel information
    fprintf(logFile, ['%.2f - Channel removed by Automatic method(PREP): '], toc);
    for ChanList = 1:size(TMP_list,2)
      fprintf(logFile, '%s(%d)\t', RMchan2{ChanList},REJchan_2(ChanList)); %  \t make space
    end
       fprintf(logFile, '\n');     
else
    fprintf(logFile, ['%.2f - Channel removed by Automatic method(PREP): None \n'], toc);
    REJchan_2=[];
end
% save automatically rejected channels information 
save('REJchan_2', 'REJchan_2');

% save current EEG temporally for compare channels between ch_rej2 and ch_rej3
EEGdummyIn = EEG;
clear bad_colors chanlist ChanList EEGdummyOut referenceIn referenceOut TMP_list indelec1 indelec2

%% Automatic Continuous Rejection
% eeglab - Tools - Automatic Continuous Rejection
% if you have trouble that 'Automatic Continuous Rejection' result in 
% data without events(triggers) except 'boundary',
% edit pop_rejcont.m file as follow.... in about line 312~313
%        else
%             NEWEEG = pop_select(EEG, 'nopoint', round(selectedregions));
%             EEG = NEWEEG;
%         end;
%     else
%         EEG = [];
%     end;

% open 'pop_rejcont' and edit line 294 to adjust 'winlength' (50 -> 30)
EEG = pop_rejcont(EEG, 'elecrange',[1:size(EEG.data,1)-3] ,'freqlimit',[10 100] ,'threshold',12,'epochlength',0.5,'contiguous',4,'addlength',0.2,'taper','hamming','eegplot','on');

% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('Continuous Rejection', 'Confirm');

% set name current file
EEG.setname = [pathsplit{1,end} '_REJcont2'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,end} '_REJcont2'];
eeglab redraw;

% Save rejected continuous data by automatic rejection
REJcont_2 = []; 
if ~isempty(TMPREJ)
REJcont_2 = TMPREJ(:,1:2); end
save ('REJcont_2','REJcont_2');
fprintf(logFile, ['%.2f - Continuous data rejection by Automatic rejection: REJcont_2.mat \n'], toc);
TMPREJ = [];

%% Find Bad Channels second time
% find bad channels by using "robust reference" in PREP
EEGdummyIn = EEG; 
referenceIn = struct('referenceChannels',[1:size(EEG.data,1)-3],'evaluationChannels',[1:size(EEG.data,1)-3],...
    'rereference',[1:size(EEG.data,1)-3],'referenceType','robust','channelLocations',EEG.chanlocs,...
    'channelInfo',EEG.chaninfo,'srate',EEG.srate,...
    'robustDeviationThreshold',7.5,'highFrequencyNoiseThreshold',4.5,'correlationThreshold',0.3,'ransacOff',true);
[EEGdummyOut, referenceOut] = performReference(EEGdummyIn,referenceIn);
disp('   PREP bad channels search is done!');

% ----------------------------------------------
% additional find High-frequency noisy channel
% find channel number F7, F5, or F3
TMP_start = find(strcmp({EEG.chanlocs.labels},'F7')==1);
if isempty(TMP_start), TMP_start = find(strcmp({EEG.chanlocs.labels},'F5')==1);
if isempty(TMP_start), TMP_start = find(strcmp({EEG.chanlocs.labels},'F3')==1); 
if isempty(TMP_start), TMP_start = find(strcmp({EEG.chanlocs.labels},'F1')==1); 
end, end, end

% find channel number P8, P6, or P4
TMP_end = find(strcmp({EEG.chanlocs.labels},'TP8')==1);
if isempty(TMP_end), TMP_end = find(strcmp({EEG.chanlocs.labels},'CP6')==1);
if isempty(TMP_end), TMP_end = find(strcmp({EEG.chanlocs.labels},'CP4')==1);    
if isempty(TMP_end), TMP_end = find(strcmp({EEG.chanlocs.labels},'CP2')==1);    
end, end, end

% find bad channels additionaly using 'Automatic channel rejection'
TMP_list=[];
[TMP_EEG, indelec1]= pop_rejchan(EEG, 'elec',[1:EEG.nbchan-3] ,'threshold',3,'norm','on','measure','spec','freqrange',[0 1] );
[TMP_EEG, indelec2]= pop_rejchan(EEG, 'elec',[TMP_start:TMP_end] ,'threshold',2.5,'norm','on','measure','spec','freqrange',[50 100] );
indelec2 = indelec2 + TMP_start - 1; %adjust channel number because it start in TMP_start(F7)
TMP_list = unique([TMP_list referenceOut.badChannels.all indelec1 indelec2]);
% find removed channels name
[~, RMchan3] = eeg_decodechan(EEG.chanlocs, TMP_list);
clear TMP_EEG TMP_start TMP_end 
disp(' =============================== ');
disp('Find Bad Channels process is done! ');
disp(' =============================== ');

% plotting EEG signal panel to see 'reject' channels
if ~isempty(TMP_list)
    chanlist = [1:EEG.nbchan-3];
    bad_colors = cell(1,length(chanlist)); bad_colors(:) = { 'k' };
    bad_colors(TMP_list) = { 'r' }; bad_colors = bad_colors(end:-1:1);

    tmpcom = [ 'EEGTMP = pop_select(EEG, ''nochannel'', [' num2str(chanlist(TMP_list)) ']);' ];
    tmpcom = [ tmpcom ...
                '[ALLEEG EEG CURRENTSET tmpcom] = pop_newset(ALLEEG, EEGTMP, CURRENTSET);' ...
                '   if ~isempty(tmpcom),' ... 
                '     EEG = eegh(LASTCOM, EEG);' ...
                '     eegh(tmpcom);' ...
                '     eeglab(''redraw'');' ...
                '  end; clear EEGTMP tmpcom;' ];

    eegplot(EEG.data(chanlist,:,:), 'srate', EEG.srate, 'title', 'Automatically selected channels to reject', ...
                 'limits', [EEG.xmin EEG.xmax]*1000, 'color', bad_colors(end:-1:1), 'eloc_file', EEG.chanlocs(chanlist), 'command', tmpcom,...
                 'winlength',30,'spacing',40);
  
    % make 'Continue button' to wait until channel rejection process is done
    SPREP2020_Nextbutton('2nd Finding Bad channels', 'Confirm');
    
    % set name current file
    EEG.setname = [name '_REJchan3'];
    ALLEEG(CURRENTSET).setname = [name '_REJchan3'];
    eeglab redraw;
    
    % find removed channels number in original data
    for j = 1:size(RMchan3,2)
    REJchan_3(j) = find(strcmp({KeepEEG.chanlocs.labels}, RMchan3{j})==1); end
    % log file automatically removed channel information
    fprintf(logFile, ['%.2f - Channel removed by Automatic method(PREP): '], toc);
    for ChanList = 1:size(TMP_list,2)
      fprintf(logFile, '%s(%d)\t', RMchan3{ChanList},REJchan_3(ChanList)); %  \t make space
    end
       fprintf(logFile, '\n');     
else
    fprintf(logFile, ['%.2f - Channel removed by Automatic method(PREP): None \n'], toc);
    REJchan_3=[];
end
% save automatically rejected channels information 
save('REJchan_3', 'REJchan_3');

% save current EEG temporally for compare channels between ch_rej2 and ch_rej3
EEGdummyIn = EEG;
clear bad_colors chanlist ChanList EEGdummyOut referenceIn referenceOut TMP_list indelec1 indelec2

%% Continuous data rejection second time
% open 'pop_rejcont' and edit line 294 to adjust 'winlength' (50 -> 30)
EEG = pop_rejcont(EEG, 'elecrange',[1:size(EEG.data,1)-3] ,'freqlimit',[10 100] ,'threshold',12,'epochlength',0.5,'contiguous',4,'addlength',0.2,'taper','hamming','eegplot','on');

% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('2nd Continuous Rejection', 'Confirm');

% set name current file
EEG.setname = [pathsplit{1,end} '_REJcont3'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,end} '_REJcont3'];
eeglab redraw;

% Save rejected continuous data by automatic rejection
REJcont_3 = []; 
if ~isempty(TMPREJ)
REJcont_3 = TMPREJ(:,1:2); end
save ('REJcont_3','REJcont_3');
fprintf(logFile, ['%.2f - Continuous data rejection by last selection: REJcont_3.mat \n'], toc);
TMPREJ = [];

%% Additional channel rejection (manually)
chanlist = [1:EEG.nbchan-3];
eegplot(EEG.data(chanlist,:,:), 'srate', EEG.srate, 'title', 'Select Additional channels to reject', ...
                 'eloc_file', EEG.chanlocs(chanlist),'winlength',30,'spacing',40);

EEG = pop_select(EEG);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 

% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('3rd Finding Bad Channels', 'Confirm');

% set name current file
EEG.setname = [name '_REJchan4'];
ALLEEG(CURRENTSET).setname = [name '_REJchan4'];
eeglab redraw;


% Save 'Remove Extra bad channels(if exist)'
RMchan4 = []; CHAN_whole = {}; CHAN_remove = {}; REJchan_4 = [];

if EEGdummyIn.nbchan ~= EEG.nbchan
    % find removed channel information
    for j = 1:size(EEGdummyIn.chanlocs,2)
    CHAN_whole{1,j} = EEGdummyIn.chanlocs(j).labels; end
    for k = 1:size(EEG.chanlocs,2)
    CHAN_remove{1,k} = EEG.chanlocs(k).labels; end
    %find channel name
    RMchan4 = setdiff(CHAN_whole, CHAN_remove,'stable'); 
    for j = 1:size(RMchan4,2)
    REJchan_4(j) = find(strcmp({KeepEEG.chanlocs.labels}, RMchan4{j})==1); end
    % there is a function, eeg_decodechan. remain the code for reference.
    % [~, RMChanlist] = eeg_decodechan(EEGdummyIn.chanlocs, indelec1);
    fprintf(logFile, ['%.2f - Channel removed by manual selection(1): '], toc);
    for ChanList = 1:size(REJchan_4,2)
    fprintf(logFile, '%s(%d)\t', RMchan4{ChanList},REJchan_4(ChanList)); %  \t make space
    end
    fprintf(logFile, '\n');
else
    fprintf(logFile, ['%.2f - Channel removed by manual selection(1): None \n'], toc);
    REJchan_4=[];
end
    % save channel list to mat file
    save('REJchan_4', 'REJchan_4');
 
clear j k CHAN_whole CHAN_remove
clear REJchan_1 REJchan_2 REJchan_3 REJchan_4 REJcont_1 REJcont_2 REJcont_3

%% Save EEG dataset
cd('..') % go to 'up'folder
savename = [name '_reject.set']
EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% log save process
fprintf(logFile, ['%.2f - Save EEG dataset(rejection done): %s \n'], toc, savename);
clear savename

%% Interpolation channels
Chan_nb = EEG.nbchan;
EEG = pop_interp(EEG, KeepEEG.chanlocs, 'spherical');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [name '_interpol'];
ALLEEG(CURRENTSET).setname = [name '_interpol'];
EEG = eeg_checkset( EEG ); eeglab redraw;
% log
Chan_nb_int = EEG.nbchan - Chan_nb;
fprintf(logFile, ['%.2f - Channels are interporated: (%d)n \n'], toc, Chan_nb_int);
clear Chan_nb Chan_nb_int

%% Keeping and remove VEO, HEO, EKG channel data temporally
Data_HEO = EEG.data(EEG.nbchan-2,:);
Data_VEO = EEG.data(EEG.nbchan-1,:);
Data_EKG = EEG.data(EEG.nbchan,:);

% remove EKG channel
EEG = pop_select( EEG,'nochannel',{'HEO','VEO','EKG'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [name '_withoutEKG'];
ALLEEG(CURRENTSET).setname = [name '_withoutEKG'];
EEG = eeg_checkset( EEG ); eeglab redraw;
% log
fprintf(logFile, ['%.2f - HEO(%d),VEO(%d),EKG(%d) channel is removed temporally \n'], toc, EEG.nbchan-2, EEG.nbchan-1, EEG.nbchan);

%% Run AMICA using calculated data rank with 'pcakeep' option
% Note : rank function built-in in MATLAB doesn't consider interpolated
% channels. So, use calculation of datarank follows.

% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('AMICA', 'Run AMICA');

% calculate dataRank 
% dataRank = rank(double(EEG.data'));  <- this function sometimes not
% working properly!
covarianceMatrix = cov(double(EEG.data(1:EEG.nbchan,1:3000))');
[E, D] = eig (covarianceMatrix);
rankTolerance = 1e-7;
dataRank = sum (diag (D) > rankTolerance)
clear covarianceMatrix E D rankTolerance

% log
fprintf(logFile, ['%.2f - ICA calculation with Rank : %d \n'], toc, dataRank);

% EEG = pop_runica(EEG, 'extended',1,'pca',dataRank,'interupt','on');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% run AMICA
% input EEG.data(1:end-1,:) for except EKG channel
dataName = 'AMICA'; 
runamica15(EEG.data(1:end,:), 'num_chans', EEG.nbchan,...
    'outdir', [pwd '\' dataName],'pcakeep', dataRank, 'num_models', 1,...
    'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1,'max_threads',p.AMICAthreads);
%if you not write fixed folder location, ex)['\amicaResults\',dataName]
%program select current drive(ex. C:\) -> C:\amicaResults\dataName

eeglab redraw;

fprintf(logFile, ['%.2f - ICA calculation is done \n'], toc);

%% load AMICA file
% this process require AMICA folder calculated already
dataName = 'AMICA'; 
EEG.etc.amica  = loadmodout15([pwd '\' dataName]);
EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
EEG.icaweights = EEG.etc.amica.W;
EEG.icasphere  = EEG.etc.amica.S;
EEG = eeg_checkset(EEG, 'ica');
EEG.setname = [pathsplit{1,end} '_AMICA'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,end} '_AMICA'];
eeglab redraw;

% log
fprintf(logFile, ['%.2f - AMICA load \n'], toc);
disp(' AMICA Process is done! ');
clear dataName 

%% Restore EKG channel
% HEO
EEG=pop_chanedit(EEG, 'append',EEG.nbchan,'changefield',{EEG.nbchan+1 'labels' 'HEO'});
EEG.nbchan = EEG.nbchan+1;
EEG.chanlocs(EEG.nbchan).labels = 'HEO'; EEG.chanlocs(EEG.nbchan).urchan = 63;
EEG.chaninfo.nodatchans = [];
EEG.data(EEG.nbchan,:) = Data_HEO;

% VEO
EEG=pop_chanedit(EEG, 'append',EEG.nbchan,'changefield',{EEG.nbchan+1 'labels' 'VEO'});
EEG.nbchan = EEG.nbchan+1;
EEG.chanlocs(EEG.nbchan).labels = 'VEO'; EEG.chanlocs(EEG.nbchan).urchan = 64;
EEG.chaninfo.nodatchans = [];
EEG.data(EEG.nbchan,:) = Data_VEO;

% EKG
EEG=pop_chanedit(EEG, 'append',EEG.nbchan,'changefield',{EEG.nbchan+1 'labels' 'EKG'});
EEG.nbchan = EEG.nbchan+1;
EEG.chanlocs(EEG.nbchan).labels = 'EKG'; EEG.chanlocs(EEG.nbchan).urchan = 65;
EEG.chaninfo.nodatchans = [];
EEG.data(EEG.nbchan,:) = Data_EKG;

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [name '_withEKG'];
ALLEEG(CURRENTSET).setname = [name '_withEKG'];
EEG = eeg_checkset( EEG ); eeglab redraw;

% log
fprintf(logFile, ['%.2f - HEO(%d),VEO(%d),EKG(%d) channel is restored \n'], toc, EEG.nbchan-2, EEG.nbchan-1, EEG.nbchan);

%% Save AMICA set file
% save AMICA set file
if p.AMICAsave == 1
    savename = [name '_AMICA.set']
    EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % log save process
    fprintf(logFile, ['%.2f - Save EEG dataset(AMICA): %s \n'], toc, savename);
    clear savename
end

%% Epoching

EEG = pop_epoch( EEG, p.Epochtrigger, p.TMP_epochrange, 'newname', [name '_Epoch'], 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG = eeg_checkset( EEG );
eeglab redraw;
% log save process
fprintf(logFile, ['%.2f - Epoching(long) temporally : %.1fs ~ %.1fs \n'], toc, p.TMP_epochrange(1),p.TMP_epochrange(2));
%% SASICA 
% EOG = auto 4, [UVEO LVEO], [RHEO LHEO]
% ECG = auto 4, [ECG]
% Autocorrelation = auto 2
% Focal components = auto 3
% Focal trial activity = auto 3

% When SASICA error like...
% SASICA>struct2vararg (line 759) input should be 1x1 structure array
% command follow code and select any dataset to reset SASICA parameters
%  >> SASICA resetprefs
% or reset preference setting in hidden MATLAB setting folder
% ex)C:\Users\admin\AppData\Roaming\MathWorks\MATLAB\R2016a\matlabprefs.mat
% if you open the file, you can find SASICA structure, remove it and
% save that the mat file.

% Or, follow code insert in SASICA.m file
% at line before lines using ' struct2vararg(cfg) '
% that, insert following code in about line 176 (ver 1.3.4)(after getpref function
% and end function. In other words, before struct2vararg function )
% In ver 1.3.7, (put the following code about line 630 (before making
% varialbe 'com')
% check cfg structure field to clear out unnecessary field
% --------------------------------------------------------------
% deletefield = {'setname','filename','filepath','subject','group',...
%     'condition','session','comments','nbchan','trials','pnts',...
%     'srate','xmin','xmax','times','data','icaact','icawinv',...
%     'icasphere','icaweights','icachansind','chanlocs','urchanlocs',...
%     'chaninfo','ref','event','urevent','eventdescription','epoch',...
%     'epochdescription','reject','stats','specdata','specicaact',...
%     'splinefile','icasplinefile','dipfit','history','saved',...
%     'etc','datfile'};
% hasfield = isfield(cfg,deletefield);
% if hasfield
%     cfg = rmfield(cfg,deletefield(find(hasfield)));
% else
%     disp('Stored preference setting has no problem');
% end
% clear deletefield hasfield;
% --------------------------------------------------------------
SASICA;

% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('After SASICA', 'Execute ICLable');

% EEGLAB ICLable function
EEG = pop_iclabel(EEG);

% Plotting ICLables
pop_viewprops(EEG, 0);
pop_selectcomps(EEG, [1:size(EEG.icaweights,1)] );

% make 'Continue button' to wait until channel rejection process is done
SPREP2020_Nextbutton('Remove ICA', 'Remove ICA');

% find selected components
REJ_comp = find(EEG.reject.gcompreject == true);
% back to continuous data
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',CURRENTSET-1,'study',0); 
EEG = eeg_checkset( EEG ); eeglab redraw;

% remove selected components
EEG = pop_subcomp( EEG, REJ_comp, 1); % selected component numbers are stored in EEG.reject.gcompreject
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',[name '_ICA remove'],'gui','off'); 
eeglab redraw;

% log save process
fprintf(logFile, ['%.2f - SASICA & ICLable - removed components: '], toc);
for j = 1:size(REJ_comp,2)
    fprintf(logFile, '%d  ', REJ_comp(j)); 
end
fprintf(logFile, '\n');
clear REJ_comp j


%% Save eeglab file
savename = [name '_conti_Final.set'];
EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
fprintf(logFile, ['%.2f - Save EEG dataset(Continuous Final): %s \n\n'], toc, savename);

fclose(logFile); %stop logging
disp(' =============================== ');
disp('    All preprocessing id done! ');
disp(' =============================== ');

