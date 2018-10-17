% Preprocessing file after two CIT files were merged
% made by SIU

%% Get name and subject number
% divide current EEG file name for search physiology file which has same
% name
[pathstr, name, ext] = fileparts(EEG.comments);
% divide path(folder)
pathsplit = strsplit(pathstr, '\');
% current folder location
Rootpath = pwd;
clear pathstr ext
%% Keeping original EEG data
% (for interpolation later)
EEG.setname = [pathsplit{1,4} '_CIT_merge']; eeglab redraw;
KeepEEG = EEG; 
tic;
%% Start Log
cd(pathsplit{1,4}); mkdir('Log'); cd('Log');
logFile = fopen([pathsplit{1,4} '_Preprocess_log.txt'], 'at+');
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
EEG.setname = [pathsplit{1,4} '_CIT_merge_REJchan1'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJchan1'];

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
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

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
% set name current file
EEG.setname = [pathsplit{1,4} '_CIT_merge_REJcont1'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJcont1'];
eeglab redraw;

% Save rejected continuous data by visual inspection
REJcont_1 = [];
if ~isempty(TMPREJ)
REJcont_1 = TMPREJ(:,1:2); end
save ('REJcont_1','REJcont_1');
fprintf(logFile, ['%.2f - Continuous data rejection by Visual inspection: REJcont_1.mat \n'], toc);
TMPREJ = [];

%% Automatic Continuous Rejection
% eeglab - Tools - Automatic Continuous Rejection
% if you have trouble that 'Automatic Continuous Rejection' result in 
% data without events(triggers) except 'boundary',
% edit eeg_rejcont.m file as follow.... in about line 298~305
%        else
%             NEWEEG = pop_select(EEG, 'nopoint', round(selectedregions));
%             EEG = NEWEEG;
%         end;
%     else
%         EEG = [];
%     end;

% open 'pop_rejcont' and edit line 294 to adjust 'winlength' (50 -> 30)
EEG = pop_rejcont(EEG, 'elecrange',[1:size(EEG.data,1)-5] ,'freqlimit',[10 100] ,'threshold',10,'epochlength',0.5,'contiguous',4,'addlength',0.2,'taper','hamming','eegplot','on');

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
% set name current file
EEG.setname = [pathsplit{1,4} '_CIT_merge_REJcont2'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJcont2'];
eeglab redraw;

% Save rejected continuous data by automatic rejection
REJcont_2 = []; 
if ~isempty(TMPREJ)
REJcont_2 = TMPREJ(:,1:2); end
save ('REJcont_2','REJcont_2');
fprintf(logFile, ['%.2f - Continuous data rejection by Automatic rejection: REJcont_2.mat \n'], toc);
TMPREJ = [];

%% Find Bad Channels
% find bad channels by using "robust reference" in PREP
EEGdummyIn = EEG; 
referenceIn = struct('referenceChannels',[1:size(EEG.data,1)-5],'evaluationChannels',[1:size(EEG.data,1)-5],...
    'rereference',[1:size(EEG.data,1)-5],'referenceType','robust','channelLocations',EEG.chanlocs,...
    'channelInfo',EEG.chaninfo,'srate',EEG.srate,...
    'robustDeviationThreshold',7,'highFrequencyNoiseThreshold',3.5,'correlationThreshold',0.3,'ransacOff',true);
[EEGdummyOut, referenceOut] = performReference(EEGdummyIn,referenceIn);
disp('   PREP bad channels search is done!');

% ----------------------------------------------
% additional find High-frequency noisy channel
% find channel number F7, F5, or F3
TMP_start = find(strcmp({EEG.chanlocs.labels},'F7')==1);
if isempty(TMP_start)
    TMP_start = find(strcmp({EEG.chanlocs.labels},'F5')==1);
if isempty (TMP_start) 
    TMP_start = find(strcmp({EEG.chanlocs.labels},'F3')==1); 
if isempty (TMP_start) 
TMP_start = find(strcmp({EEG.chanlocs.labels},'F1')==1); 
end, end, end

% find channel number P8, P6, or P4
TMP_end = find(strcmp({EEG.chanlocs.labels},'TP8')==1);
if isempty(TMP_end)
    TMP_end = find(strcmp({EEG.chanlocs.labels},'CP6')==1);
if isempty (TMP_end)
    TMP_end = find(strcmp({EEG.chanlocs.labels},'CP4')==1);    
if isempty (TMP_end)
    TMP_end = find(strcmp({EEG.chanlocs.labels},'CP2')==1);    
end, end, end

% find bad channels additionaly using 'Automatic channel rejection'
TMP_list=[];
[TMP_EEG, indelec1]= pop_rejchan(EEG, 'elec',[1:EEG.nbchan-5] ,'threshold',3.5,'norm','on','measure','spec','freqrange',[0 1] );
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
    chanlist = [1:EEG.nbchan-5];
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
    F_check = figure;
    set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
    B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
    uiwait(F_check); close(F_check);    
    clear B_conti F_check
    % set name current file
    EEG.setname = [pathsplit{1,4} '_CIT_merge_REJchan2'];
    ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJchan2'];
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

%% Additional channel rejection (manually)
chanlist = [1:EEG.nbchan-5];
eegplot(EEG.data(chanlist,:,:), 'srate', EEG.srate, 'title', 'Select Additional channels to reject', ...
                 'eloc_file', EEG.chanlocs(chanlist),'winlength',30,'spacing',40);

EEG = pop_select(EEG);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
% set name current file
EEG.setname = [pathsplit{1,4} '_CIT_merge_REJchan3'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJchan3'];
eeglab redraw;


% Save 'Remove Extra bad channels(if exist)'
RMchan3 = []; CHAN_whole = {}; CHAN_remove = {}; REJchan_3 = [];

if EEGdummyIn.nbchan ~= EEG.nbchan
    % find removed channel information
    for j = 1:size(EEGdummyIn.chanlocs,2)
    CHAN_whole{1,j} = EEGdummyIn.chanlocs(j).labels; end
    for k = 1:size(EEG.chanlocs,2)
    CHAN_remove{1,k} = EEG.chanlocs(k).labels; end
    %find channel name
    RMchan3 = setdiff(CHAN_whole, CHAN_remove,'stable'); 
    for j = 1:size(RMchan3,2)
    REJchan_3(j) = find(strcmp({KeepEEG.chanlocs.labels}, RMchan3{j})==1); end
    % there is a function, eeg_decodechan. remain the code for reference.
    % [~, RMChanlist] = eeg_decodechan(EEGdummyIn.chanlocs, indelec1);
    fprintf(logFile, ['%.2f - Channel removed by manual selection(1): '], toc);
    for ChanList = 1:size(REJchan_3,2)
    fprintf(logFile, '%s(%d)\t', RMchan3{ChanList},REJchan_3(ChanList)); %  \t make space
    end
    fprintf(logFile, '\n');
else
    fprintf(logFile, ['%.2f - Channel removed by manual selection(1): None \n'], toc);
    REJchan_3=[];
end
    % save channel list to mat file
    save('REJchan_3', 'REJchan_3');

    
%% last continuous data rejection
% open 'pop_rejcont' and edit line 294 to adjust 'winlength' (50 -> 30)
EEG = pop_rejcont(EEG, 'elecrange',[1:size(EEG.data,1)-5] ,'freqlimit',[10 100] ,'threshold',12,'epochlength',0.5,'contiguous',4,'addlength',0.2,'taper','hamming','eegplot','on');

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
% set name current file
EEG.setname = [pathsplit{1,4} '_CIT_merge_REJcont3'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJcont3'];
eeglab redraw;

% Save rejected continuous data by automatic rejection
REJcont_3 = []; 
if ~isempty(TMPREJ)
REJcont_3 = TMPREJ(:,1:2); end
save ('REJcont_3','REJcont_3');
fprintf(logFile, ['%.2f - Continuous data rejection by last selection: REJcont_3.mat \n'], toc);
TMPREJ = [];

%% Additional channel rejection (manually) one more!
% save current EEG temporally for compare channels between ch_rej4
EEGdummyIn = EEG;

chanlist = [1:EEG.nbchan-5];
eegplot(EEG.data(chanlist,:,:), 'srate', EEG.srate, 'title', 'Select Additional channels to reject', ...
                 'eloc_file', EEG.chanlocs(chanlist),'winlength',30,'spacing',40);

EEG = pop_select(EEG);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
% set name current file
EEG.setname = [pathsplit{1,4} '_CIT_merge_REJchan4'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_REJchan4'];
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
    fprintf(logFile, ['%.2f - Channel removed by manual selection(2): '], toc);
    for ChanList = 1:size(REJchan_4,2)
    fprintf(logFile, '%s(%d)\t', RMchan4{ChanList},REJchan_4(ChanList)); %  \t make space
    end
    fprintf(logFile, '\n');
else
    fprintf(logFile, ['%.2f - Channel removed by manual selection(2): None \n'], toc);
    REJchan_4=[];
end
    % save channel list to mat file
    save('REJchan_4', 'REJchan_4');
    
clear j k CHAN_whole CHAN_remove
clear REJchan_1 REJchan_2 REJchan_3 REJchan_4 REJcont_1 REJcont_2 REJcont_3
%% Save EEG dataset
cd('..') % go to 'up'folder
savename = [pathsplit{1,4} '_CIT_[0_5Hz]Reject.set']
EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% log save process
fprintf(logFile, ['%.2f - Save EEG dataset(rejection done): %s \n'], toc, savename);
clear savename

%% Interpolation channels
Chan_nb = EEG.nbchan;
EEG = pop_interp(EEG, KeepEEG.chanlocs, 'spherical');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [pathsplit{1,4} '_CIT_merge_interpol'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_interpol'];
EEG = eeg_checkset( EEG ); eeglab redraw;
% log
Chan_nb_int = EEG.nbchan - Chan_nb;
fprintf(logFile, ['%.2f - Channels are interporated: (%d)n \n'], toc, Chan_nb_int);
clear Chan_nb Chan_nb_int

%% Keeping and remove ECG channel data temporally
Data_ECG = EEG.data(EEG.nbchan,:);
% remove ECG channel
EEG = pop_select( EEG,'nochannel',{'ECG'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [pathsplit{1,4} '_CIT_merge_withoutECG'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_withoutECG'];
EEG = eeg_checkset( EEG ); eeglab redraw;
% log
fprintf(logFile, ['%.2f - ECG channel(%d) is removed temporally \n'], toc, EEG.nbchan);

%% Run AMICA using calculated data rank with 'pcakeep' option
% Sometimes, AMICA ignore dataRank.. I guess EEGLAB pop_newset function is
% late for next step, so CURRENTSET - 1 goes into AMICA procedure
% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Start AMICA','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','AMICA(ICA)',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check

% calculate dataRank 
% dataRank = rank(double(EEG.data'));  <- this function sometimes not
% working properly!
covarianceMatrix = cov(double(EEG.data(1:EEG.nbchan,1:3000))');
[E, D] = eig (covarianceMatrix);
rankTolerance = 1e-7;
dataRank = sum (diag (D) > rankTolerance);
clear covarianceMatrix E D rankTolerance

% log
fprintf(logFile, ['%.2f - AMICA calculation with Rank : %d \n'], toc, dataRank);

% run AMICA
% input EEG.data(1:end-1,:) for except ECG channel
dataName = 'AMICA'; 
runamica15(EEG.data(1:end,:), 'num_chans', EEG.nbchan,...
    'outdir', [pwd '\' dataName],'pcakeep', dataRank, 'num_models', 1,...
    'do_reject', 1, 'numrej', 15, 'rejsig', 3, 'rejint', 1,'max_threads',4);
%if you not write fixed folder location, ex)['\amicaResults\',dataName]
%program select current drive(ex. C:\) -> C:\amicaResults\dataName

disp(' =============================== ');
disp('    All preprocessing id done! ');
disp(' =============================== ');

fprintf(logFile, ['%.2f - AMICA calculation is done \n'], toc);
fclose(logFile); %stop logging