% Preprocessing Semi-automatic for 0.1hz filtered data
% made by SIU

%% load data .cnt file
% In here, you can load .cnt file in eeglab GUI.

% eeglab;
% setting option not to use single precision.('option_single')
% pop_editoptions('option_storedisk', 0, 'option_savetwofiles', 1, 'option_saveversion6', 1,...
%     'option_single', 0, 'option_memmapdata', 0, 'option_eegobject', 0, 'option_computeica', 1,...
%     'option_scaleicarms', 1, 'option_rememberfolder', 1, 'option_donotusetoolboxes', 0,...
%     'option_checkversion', 1, 'option_chat', 0);
%% set parameter
Range_longepoch         = [-1 1.7];         % second
Range_microdetect       = [-200 1200];      % ms
p.Usemicro              = 0;                % if 1 = Doing microsaccade procedure // 0 = not (If your target is P3, I recommend 0 )
p.Usemicrorange         = 0;                % if 1 = using defined microdetect range // 0 = whole epoch range

p.AMICAsave             = 0;                % if 1 = save .set file // 0 = no save

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
EEG.setname = [pathsplit{1,4} '_CIT_merge0.1']; eeglab redraw;
KeepEEG = EEG; EEGdummyIn = EEG; 
tic;
%% Start log
tic;
cd(pathsplit{1,4}); mkdir('Log'); cd('Log');
logFile = fopen([pathsplit{1,4} '_Preprocess_log.txt'], 'at+');
fprintf(logFile, ['\n===================================== \n']);
fprintf(logFile, [' START: ' datestr(now) '\n']);
fprintf(logFile, [' FILE : ' name ' \n']);
fprintf(logFile, [' ¡Ø This file is filter at 0.1 and load AMICA \n \n']);

%% Load Rejected channel info & Rejected continuous data information
load REJchan_1.mat; load REJchan_2.mat; load REJchan_3.mat;  load REJchan_4.mat;
load REJcont_1.mat; load REJcont_2.mat; load REJcont_3.mat;
REJ_list = sort([REJchan_1, REJchan_2, REJchan_3 REJchan_4]);
fprintf(logFile, ['%.2f - Load Rejection info (Channel, Time points \n'], toc);
%% Reject channels
% find channel names
if ~isempty(REJ_list)
    [~, RMChanlist] = eeg_decodechan(KeepEEG.chanlocs, REJ_list);
    EEG = pop_select( EEG,'nochannel',REJ_list);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
    % log reject channel
    fprintf(logFile, ['%.2f - Channel rejection : '], toc);
    for ChanList = 1:size(RMChanlist,2)
    fprintf(logFile, '%s(%d)\t', RMChanlist{ChanList},REJ_list(ChanList)); %  \t make space
    end
    fprintf(logFile, '\n');
else
    fprintf(logFile, ['%.2f - Channel rejection : None \n'], toc);
end

clear REJ_list RMChanlist ChanList REJchan_1 REJchan_2 REJchan_3 REJchan_4
%% Reject continuous segments
% first sements rejection
EEG = pop_select(EEG, 'nopoint', REJcont_1);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
fprintf(logFile, ['%.2f - Continuous data segments rejection : First \n'], toc);
% second sements rejection
EEG = pop_select(EEG, 'nopoint', REJcont_2);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
fprintf(logFile, ['%.2f - Continuous data segments rejection : Second \n'], toc);
% third sements rejection
EEG = pop_select(EEG, 'nopoint', REJcont_3);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
fprintf(logFile, ['%.2f - Continuous data segments rejection : Third \n'], toc);

%recheck EEGLAB
EEG = eeg_checkset( EEG );
eeglab redraw;
clear REJcont_1 REJcont_2 REJcont_3

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

%% load AMICA file
% this process require AMICA folder calculated already
cd('..');
dataName = 'AMICA'; 
EEG.etc.amica  = loadmodout15([pwd '\' dataName]);
EEG.etc.amica.S = EEG.etc.amica.S(1:EEG.etc.amica.num_pcs, :); % Weirdly, I saw size(S,1) be larger than rank. This process does not hurt anyway.
EEG.icaweights = EEG.etc.amica.W;
EEG.icasphere  = EEG.etc.amica.S;
EEG = eeg_checkset(EEG, 'ica');
EEG.setname = [pathsplit{1,4} '_CIT_merge_AMICA'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_AMICA'];
eeglab redraw;

% log
fprintf(logFile, ['%.2f - AMICA load \n'], toc);
disp(' AMICA Process is done! ');
clear dataName 
%% Restore ECG channel & save _AMICA.set
EEG=pop_chanedit(EEG, 'append',EEG.nbchan,'changefield',{EEG.nbchan+1 'labels' 'ECG'});
EEG.nbchan = EEG.nbchan+1;
EEG.chanlocs(EEG.nbchan).labels = 'ECG'; EEG.chanlocs(EEG.nbchan).urchan = 69;
EEG.chaninfo.nodatchans = [];
EEG.data(EEG.nbchan,:) = Data_ECG;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [pathsplit{1,4} '_CIT_merge_withECG'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_merge_withECG'];
EEG = eeg_checkset( EEG ); eeglab redraw;
% log
fprintf(logFile, ['%.2f - ECG channel(%d) is restored \n'], toc, EEG.nbchan);

% save AMICA set file
if p.AMICAsave == 1
    savename = [pathsplit{1,4} '_CIT_[0_1Hz]AMICA.set']
    EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    % log save process
    fprintf(logFile, ['%.2f - Save EEG dataset(AMICA): %s \n'], toc, savename);
    clear savename
end
%% Epoching
EEG = pop_epoch( EEG, {  '10'  '20'  '30'  '40'  '50'  '60'  '70'  }, Range_longepoch, 'newname', [pathsplit{1,4} '_CIT_epoch'], 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG = eeg_checkset( EEG );
eeglab redraw;
% log save process
fprintf(logFile, ['%.2f - Epoching(long) for wavelet: %.1fs ~ %.1fs \n'], toc, Range_longepoch(1),Range_longepoch(2));
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
% that, insert following code in about line 169
% check cfg structure field to clear out unnecessary field
% --------------------------------------------------------------
% deletefield = {'setname','filename','filepath','subject','group',...
%     'condition','session','comments','nbchan','trials','pnts',...
%     'srate','xmin','xmax','times','data','icaact','icawinv',...
%     'icasphere','icaweights','icachansind','chanlocs','urchanlocs',...
%     'chaninfo','ref','event','urevent','eventdescription','epoch',...
%     'epochdescription','reject','stats','specdata','specicaact',...
%     'splinefile','icasplinefile','dipfit','history','saved',...
%     'etc','datfile'}
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
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Check','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Remove components',...
          'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check

%remove selected component
REJ_comp = find(EEG.reject.gcompreject == true);
EEG = pop_subcomp( EEG, REJ_comp, 1); % selected component numbers are stored in EEG.reject.gcompreject
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname',[pathsplit{1,4} '_CIT_epoch_SASICA'],'gui','off'); 

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Remove Comp Confirm','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Remove comp confirm',...
          'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check

% log save process
fprintf(logFile, ['%.2f - SASICA - removed components: '], toc);
for j = 1:size(REJ_comp,2)
    fprintf(logFile, '%d  ', REJ_comp(j)); 
end
fprintf(logFile, '\n');
clear REJ_comp j
%% Remove ECG channel
EEG = pop_select( EEG,'nochannel',{'ECG'});
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
EEG.setname = [pathsplit{1,4} '_CIT_epoch_SASICA_noECG'];
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_epoch_SASICA_noECG'];
EEG = eeg_checkset( EEG ); eeglab redraw;
% log
fprintf(logFile, ['%.2f - ECG channel(%d) is removed \n\n'], toc, EEG.nbchan+1);

%% Microsaccade remove (Microdetect plugin)
if p.Usemicro == 1

% copy current dataset
[ALLEEG EEG CURRENTSET] = pop_copyset( ALLEEG, CURRENTSET, CURRENTSET+1); eeglab redraw;
% Keeping number of current dataset
Num_Keepdata = CURRENTSET;

% calculate rEOG
rEOGchan = [EEG.nbchan-3:EEG.nbchan];
if p.Usemicrorange == 1
    EEG = pop_rEOG(EEG, 'eyechans', rEOGchan, 'window', Range_microdetect, 'filt', 1, 'method', 1);
    fprintf(logFile, ['%.2f - Find rEOG with time range: %dms ~ %dms \n'], toc,Range_microdetect(1),Range_microdetect(2));
else
    EEG = pop_rEOG(EEG, 'eyechans', rEOGchan, 'filt', 1, 'method', 1);
    fprintf(logFile, ['%.2f - Find rEOG with time range: %dms ~ %dms \n'], toc,EEG.times(1),EEG.times(end));
end

[ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG); eeglab redraw; % this line for update ALLEEG information

% detect microsaccade  (Add saccade event check 'on')
if p.Usemicrorange == 1
    EEG = pop_detect(EEG, 'eyechans', rEOGchan, 'dataset', CURRENTSET, 'thresh', 3, 'window', Range_microdetect, 'addsacs', 1, 'normRate', 1, 'plot', 0);
    fprintf(logFile, ['%.2f - Detect microsaccades with time range: %dms ~ %dms \n'], toc,Range_microdetect(1),Range_microdetect(2));
else
    EEG = pop_detect(EEG, 'eyechans', rEOGchan, 'dataset', CURRENTSET, 'thresh', 3, 'addsacs', 1, 'normRate', 1, 'plot', 0);
    fprintf(logFile, ['%.2f - Detect microsaccades with time range: %dms ~ %dms \n'], toc,EEG.times(1),EEG.times(end));
end
eeglab redraw;

% Run ICA for detection
EEG = pop_saccICA(EEG);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'setname','Microsaccade Epochs','gui','off'); 
disp(' Microdetect ICA process is done');
fprintf(logFile, ['%.2f - Microdetect ICA process \n'], toc);

% Plot ICA component correlated with eye channels
pop_eyeCorrs(EEG, 'eyechans', rEOGchan); % command history is wrong.
pop_prop(EEG, 0); % 1 is channel, 0 is ICA components

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 270 70],'Name','Apply ICA weight','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 230 40],'String','Apply ICA weight to primary dataset',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check
REJ_micro = find(EEG.reject.gcompreject == true); %keep comp number to reject
% log
fprintf(logFile, ['%.2f - Select microsaacade candidate components: '], toc);
for j = 1:size(REJ_micro,2)
    fprintf(logFile, '%d  ', REJ_micro(j)); 
end
fprintf(logFile, '\n');

% back to original dataset, load ICA weights
Num_microdata = CURRENTSET; % keep microdetect ICA dataset number
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',Num_Keepdata,'study',0);
% load ICA weights from microdetect ICA dataset(Num_microdata)
EEG = pop_editset(EEG, 'icachansind', ['ALLEEG(' num2str(Num_microdata) ').icachansind'],...
    'icaweights', ['ALLEEG(' num2str(Num_microdata) ').icaweights'], 'icasphere', ['ALLEEG(' num2str(Num_microdata) ').icasphere']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
fprintf(logFile, ['%.2f - Load ICA weights from microsaccade dataset \n'], toc);

% pop-up components reject window
if ~isempty(REJ_micro)
    EEG.reject.gcompreject(REJ_micro) = 1; %'reject' on in current dataset
    pop_prop( EEG, 0, REJ_micro, NaN, {'freqrange' [20 100] }); % second '0' makes toggle 'rejection', otherwise no button(Nan)
    % pop_prop(EEG,0); % 1 is channel, 0 is ICA components
    [EEG REJcomp]= pop_subcomp( EEG ); % selected component numbers are stored in EEG.reject.gcompreject
else
    pop_prop( EEG, 0); % pop-up components property with null (selected)components
    [EEG REJcomp]= pop_subcomp( EEG );
end

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Continue','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Continue',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check

EEG.setname = [pathsplit{1,4} '_CIT_microdetect']; 
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_microdetect']; 
EEG = eeg_checkset( EEG ); ALLEEG = pop_delset( ALLEEG, [Num_Keepdata+1  Num_Keepdata+2] ); %remove dataset to find microsac
eeglab redraw;

% log
if ~isempty(REJcomp)
    REJcomp_tmp = strsplit(REJcomp,'['); REJcomp_tmp = strsplit(REJcomp_tmp{2}, ']');
    REJcomp_tmp = num2str(REJcomp_tmp{1}); %find rejected component number in history(REJcomp)
    fprintf(logFile, ['%.2f - Select elected microsaacade components: %s \n'], toc, REJcomp_tmp);
else
    fprintf(logFile, ['%.2f - Select elected microsaacade components: none \n'], toc);
end

% save Microdetect set file
savename = [pathsplit{1,4} '_CIT_[0_1Hz]microdetect.set']
EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
fprintf(logFile, ['%.2f - Save EEG dataset(microdetect): %s \n\n'], toc, savename);

clear REJcomp_tmp REJcomp REJ_micro Num_Keepdata Num_microdata rEOGchan j 

else
    fprintf(logFile, ['%.2f - Microdetect(microsaccade) procedure skipped \n'], toc);
end
%% Reject epochs
% change initial value : 
% open pop_rejmenu and edit line 200 as: fastif(icacomp, [ '1:' int2str(EEG.nbchan-4) ],~~~~~ 
% and line250 ~ line 277 and 

% initialization
EEG.reject.disprej = {'manual','thresh','const','jp','kurt','freq'};
TMP_zero = zeros(1, EEG.trials); 
EEG.reject.rejjp = TMP_zero; EEG.reject.rejkurt = TMP_zero; EEG.reject.rejmanual = TMP_zero;
EEG.reject.rejthresh = TMP_zero; EEG.reject.rejconst = TMP_zero; EEG.reject.rejfreq = TMP_zero;
TMP_zero = zeros(EEG.nbchan, EEG.trials);
EEG.reject.rejjpE = TMP_zero; EEG.reject.rejkurtE = TMP_zero; EEG.reject.rejmanualE = TMP_zero;
EEG.reject.rejthreshE = TMP_zero; EEG.reject.rejconstE = TMP_zero; EEG.reject.rejfreqE = TMP_zero;
clear TMP_zero

% 1-1)abnormal values : 75 // -75 // -200 // 1500 // 1:max-4 chan
EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan-4] ,-75,75,-.2,1.496,2,0);
EEG.etc.Rejepoch.thresh1 = find(EEG.reject.rejthresh);
TMP_thresh1 = EEG.reject.rejthresh; TMP_threshE1 = EEG.reject.rejthreshE;
% 1-2)abnormal values : 80 // -80 // -1200 // 2200 // 1:max-4 chan
EEG = pop_eegthresh(EEG,1,[1:EEG.nbchan-4] ,-80,80,-1.2,2.196,2,0);
EEG.etc.Rejepoch.thresh2 = find(EEG.reject.rejthresh); 
TMP_thresh2 = EEG.reject.rejthresh; TMP_threshE2 = EEG.reject.rejthreshE;
% conbine trials from two conditions
TMPt_thresh = TMP_thresh1 + TMP_thresh2; TMPt_threshE = TMP_threshE1 + TMP_threshE2;
TMPt_thresh(find(TMPt_thresh>1)) = 1; TMPt_threshE(find(TMPt_threshE>1)) = 1;
% re-put to EEG structure
EEG.reject.rejthresh = TMPt_thresh; EEG.reject.rejthreshE = TMPt_threshE;
clear TMP_thresh1 TMP_threshE1 TMP_thresh2 TMP_threshE2 TMPt_thresh TMPt_threshE

% 2)abnormal trends : 40 // 0.3 // 1:max-4 chan
EEG = pop_rejtrend(EEG,1,[1:EEG.nbchan-4] ,850,40,0.3,2,0,0); 
EEG.etc.Rejepoch.const = find(EEG.reject.rejconst);

% 3-1)abnormal spectra : 25 // -25 // 10 // 20 // 1:max-4 chan
% multiple threshold like as... : 'threshold',[-50 50;-100 25],'freqlimits',[0 2;20 40]
% but multiple threshold work as 'and' not 'or'
EEG = pop_rejspec( EEG, 1,'elecrange',[1:EEG.nbchan-4] ,'method','multitaper','threshold',[-25 25] ,'freqlimits',[10 20]); 
EEG.etc.Rejepoch.freq1 = find(EEG.reject.rejfreq);
TMP_thresh1 = EEG.reject.rejfreq; TMP_threshE1 = EEG.reject.rejfreqE;
% 3-2)abnormal spectra : 30 // -30 // 40 // 50 // 1:max-4 chan
EEG = pop_rejspec( EEG, 1,'elecrange',[1:EEG.nbchan-4] ,'method','multitaper','threshold',[-25 25] ,'freqlimits',[20 40]); 
EEG.etc.Rejepoch.freq2 = find(EEG.reject.rejfreq);
TMP_thresh2 = EEG.reject.rejfreq; TMP_threshE2 = EEG.reject.rejfreqE;
% conbine trials from two conditions
TMPt_thresh = TMP_thresh1 + TMP_thresh2; TMPt_threshE = TMP_threshE1 + TMP_threshE2;
TMPt_thresh(find(TMPt_thresh>1)) = 1; TMPt_threshE(find(TMPt_threshE>1)) = 1;
% re-put to EEG structure
EEG.reject.rejfreq = TMPt_thresh; EEG.reject.rejfreqE = TMPt_threshE;
clear TMP_thresh1 TMP_threshE1 TMP_thresh2 TMP_threshE2 TMPt_thresh TMPt_threshE

% concactenate bad epochs
EEG.etc.Rejepoch.All = horzcat(EEG.etc.Rejepoch.thresh1,EEG.etc.Rejepoch.thresh2,...
    EEG.etc.Rejepoch.const, EEG.etc.Rejepoch.freq1, EEG.etc.Rejepoch.freq2);
EEG.etc.Rejepoch.All = unique(EEG.etc.Rejepoch.All, 'sorted');
% pop_eegplot( EEG, 1, 1, 0); % EEG, 1(nonICA), 1(superpose), 1(reject button)

% Tools - Reject data epochs - Reject data (all methods)
pop_rejmenu(EEG, 1);

% make 'Continue button' to wait until channel rejection process is done
F_check = figure;
set(gcf, 'Position',[600 20 250 70],'Name','Epoch Save final file','NumberTitle','off');
B_conti = uicontrol('Position',[20 20 200 40],'String','Save final file',...
              'Callback','uiresume(gcbf)');
uiwait(F_check); close(F_check);    
clear B_conti F_check

% change setname
EEG.setname = [pathsplit{1,4} '_CIT_Rejepoch']; 
ALLEEG(CURRENTSET).setname = [pathsplit{1,4} '_CIT_Rejepoch']; 
EEG = eeg_checkset( EEG ); eeglab redraw;

% find selected 'bad epochs' info in previous(saved as 'microdetect') dataset and store.
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',CURRENTSET-1,'study',0); % back to dataset(not reject yet) temporally
FinalRejepoch.manual = find(EEG.reject.rejmanual);
FinalRejepoch.thresh = find(EEG.reject.rejthresh);
FinalRejepoch.const = find(EEG.reject.rejconst);
FinalRejepoch.freq = find(EEG.reject.rejfreq);
FinalRejepoch.All = horzcat(FinalRejepoch.manual,FinalRejepoch.thresh,...
    FinalRejepoch.const, FinalRejepoch.freq);
FinalRejepoch.All = unique(FinalRejepoch.All, 'sorted');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',CURRENTSET+1,'study',0); % come back to current dataset(Rejepoch)
% combine 'FinalRejepoch' information to current dataset(_Rejepoch)
EEG.etc.FinalRejepoch = FinalRejepoch; clear FinalRejepoch;

% log
fprintf(logFile, ['%.2f - Epoch rejection(Manual)        : %s \n'], toc, num2str(EEG.etc.FinalRejepoch.manual));
fprintf(logFile, ['%.2f - Epoch rejection(Extreme value) : %s \n'], toc, num2str(EEG.etc.FinalRejepoch.thresh));
fprintf(logFile, ['%.2f - Epoch rejection(Abnormal trend): %s \n'], toc, num2str(EEG.etc.FinalRejepoch.const));
fprintf(logFile, ['%.2f - Epoch rejection(Abnormal freq) : %s \n'], toc, num2str(EEG.etc.FinalRejepoch.freq));
fprintf(logFile, ['%.2f - Epoch rejection(ALL)           : %s \n'], toc, num2str(EEG.etc.FinalRejepoch.All));

% save reject epoch information to .mat file
cd('Log');
if ~isempty(EEG.etc.FinalRejepoch)
REJepoch = EEG.etc.FinalRejepoch; end
save ('REJepoch','REJepoch'); clear REJepoch
cd('..'); % go up to subj number folder

%% Save eeglab file
% save Microdetect set file
savename = [pathsplit{1,4} '_CIT_[0_1Hz]final.set'];
EEG = pop_saveset( EEG, 'filename',savename,'filepath',[pwd '\\']);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
fprintf(logFile, ['%.2f - Save EEG dataset(final epochs): %s \n\n'], toc, savename);

fclose(logFile); %stop logging
disp(' =============================== ');
disp('    All preprocessing id done! ');
disp(' =============================== ');