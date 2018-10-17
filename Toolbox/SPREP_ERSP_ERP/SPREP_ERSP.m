% Individual ERSP and ERP plotting and save .mat file
% made by SIU
%% initialization
addpath(pwd);
clear all;
close all;
clc;

%% set parameters
p.suffix = '_CIT_[0_1Hz]final.set';  % suffix name to load file
p.sLORETA = 0;                       % if this set as 1, save .txt file for sLORETA

p.wavelet.calcOn = 0;
p.wavelet.saveeach = 0;              % save each condition to .set file (ex._probe.set) : 1(save) 0(no save)
p.wavelet.time = [-1000 1696]        % set wavelet timepoints
p.wavelet.cycle = [2.5         10];  % set wavelet cycle number
p.wavelet.baseline = [-300 -100];    % set wavelet baseline time range
p.wavelet.freqs = [2 100];         % set wavelet frequency range (Hz range)
p.wavelet.nfreqs = 100;              % set wavelet freqs resolution
p.wavelet.erspmax = 10;              % set wavelet results dB range
p.wavelet.alpha = 0.05;              % set wavelet permutation threshold

p.ERP.saveeach = 0;                  % save each condition to .set file (ex._probe.set) : 1(save) 0(no save)
p.ERP.savemat = 1;                   % save ERP file to .mat too (for easy to access) : 1(save) 0(no)
p.ERP.lowpassfilter = 27;            % ERP epochs low-pass filter edge (for 30Hz cut-off, use 27)
p.ERP.time = [-0.2 1.5];             % ERP epochs time range (-0.2s ~ 1.5s)
p.ERP.baseline = [-200    0];        % ERP epochs baseline range
p.ERP.multiERP = 1;                  % make multiple scalp ERP picture : 1(save) 0(no)
p.ERP.topo = 1;                      % make ERP topo file : 1(save) 0(no)


% set List to auto process
Subj.num = {'3','4','5','7','8','9','10','11','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','29','30','31','32','33','34','35','36','37','38','39','41','42','43','44','45','46','47','48','49','50','51','52','54','55','56','57','58','59','60','61','62','64','65','66','67','68','69'};
Subj.group = {'G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G','G'};
Subj.target = {'지갑','지갑','지갑','지갑','지갑','지갑','지갑','지갑','지갑','지갑','목걸이','목걸이','목걸이','목걸이','목걸이','목걸이','목걸이','목걸이','목걸이','목걸이','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','만년필','지갑','지갑','지갑','지갑','지갑','지갑','지갑','목걸이','목걸이','목걸이','목걸이','목걸이','목걸이','향수','향수','향수','향수','향수','향수','향수','향수','향수','향수','향수','목걸이'};
Subj.probe = {'반지','시계','반지','반지','반지','시계','시계','시계','반지','반지','시계','반지','시계','반지','반지','시계','시계','반지','반지','시계','시계','시계','시계','반지','시계','시계','반지','시계','시계','시계','반지','시계','시계','시계','반지','시계','시계','시계','시계','반지','시계','반지','시계','시계','반지','반지','시계','반지','시계','반지','반지','시계','반지','반지','시계','시계','반지','반지','반지','시계','시계'};


for Listnum = 1:size(Subj.num,2)
% open eeglab
eeglab;
pause(2);
    
%% set code to use input
Table_code.Code = [10:10:70];
Table_code.Name = {'시계' '반지' '지갑' '목걸이' '향수' '만년필' '벨트'};

%% Get input subject number and target, probe object
% subj_num = input('subject number??   ', 's');
% subj_group = input('subject group?? (G/I)   ', 's');
% subj_target = input('Target object??   ', 's');
subj_num = Subj.num{Listnum};
subj_group = Subj.group{Listnum};
subj_target = Subj.target{Listnum};

%find target item number in Table_code
Index_target = strfind(Table_code.Name, subj_target);
Index_target = find(not(cellfun('isempty', Index_target)));
Index_tr = Table_code.Code(Index_target); 

if subj_group == 'G' || subj_group == 'g'
   subj_probe = Subj.probe{Listnum};
   %find probe item number in Table_code
   Index_probe = strfind(Table_code.Name, subj_probe);
   Index_probe = find(not(cellfun('isempty', Index_probe)));
   Index_pr = Table_code.Code(Index_probe);
else
    Index_random = randi(10);
    if Index_random > 5
        Index_pr = [10];
    else
        Index_pr = [20];
    end
end

%calculate Irrelevant item number 
Index_ir = Table_code.Code; 
for i = 1:size(Index_pr,2)
    Index_ir(find(Index_ir == Index_pr(i))) = [];
end
Index_ir(find(Index_ir == Index_tr)) = [];
Index_irrelevant = Index_ir./10;
clear Table_code_num Table_code_str i;

%% load file and Start logging
cd(subj_num); % go to subj number folder
cwd = pwd; name_load = [subj_num,p.suffix]; % name of file to load

% Start log
tic;
mkdir('Log'); cd('Log');
logFile = fopen([subj_num '_ERSP_log.txt'], 'at+');
fprintf(logFile, ['\n===================================== \n']);
fprintf(logFile, [' START: ' datestr(now) '\n']);
fprintf(logFile, [' FILE : ' name_load ' \n \n']);
cd('..'); % go to subj number folder

% load eeglab set file
EEG = pop_loadset('filename',name_load,'filepath',cwd);
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
eeglab redraw;
% log
fprintf(logFile, ['%.2f - Loading .set file: %s \n'], toc, name_load);
clear name_load

%if p.wavelet.calcOn == 1   line 138
        
%% save whole ERSP.set file which trigger changed for statistics
mkdir('TF_analysis'); cd('TF_analysis');
% re-naming trigger events
EEG = pop_selectevent( EEG, 'type',Index_tr,'renametype','Target','deleteevents','off','deleteepochs','off','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','10_CIT_ERSP_rename','gui','off'); EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_pr,'renametype','Probe','deleteevents','off','deleteepochs','off','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off'); EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_ir ,'renametype','Irrelevant','deleteevents','off','deleteepochs','off','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'overwrite','on','gui','off'); EEG = eeg_checkset( EEG );
eeglab redraw;
% log
Item_tr = Table_code.Name{Index_target}; Item_pr = Table_code.Name{Index_probe}; Item_ir = [];
for i = 1:size(Index_irrelevant,2), Item_ir = [Item_ir ' ' Table_code.Name{Index_irrelevant(i)}]; end
fprintf(logFile, ['%.2f - Change trigger number: Target( %s ) / Probe( %s ) / Irrelevant(%s) \n'], toc, Item_tr, Item_pr, Item_ir);
clear Index_target Index_probe Index_irrelevant Item_tr Item_pr Item_ir;

%save
cwd = pwd; name_save = [subj_num,'_CIT_ERSP.set'];
EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
cd('..');

ALLEEG = pop_delset( ALLEEG, [2] ); % delete renamed dataset (#2 dataset)
eeglab redraw;
% log
fprintf(logFile, ['%.2f - Save ERSP .set file: %s \n'], toc, name_save);

if p.wavelet.calcOn == 1
%% divide a dataset following triggers
EEG = pop_selectevent( EEG, 'type',Index_tr,'deleteevents','on','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Target','gui','off'); Trial_tr = EEG.trials;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_pr,'deleteevents','on','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Probe','gui','off'); Trial_pr = EEG.trials;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',1,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_ir,'deleteevents','on','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Irrelevant','gui','off'); Trial_ir = EEG.trials;
eeglab redraw;
%log
fprintf(logFile, ['%.2f - Divide ERSP file to each condition: Target(%d TR) / Probe(%d TR) / Irrelevant(%d TR) \n'], toc, Trial_tr, Trial_pr, Trial_ir);
clear Trial_tr Trial_pr Trial ir;
%% sLORETA txt file 
if p.sLORETA == 1
    mkdir('sLORETA'); cd('sLORETA');
    mkdir('ERSP'); cd('ERSP');
    % open Dataset #2 (Target) -----------------------------------------------
    mkdir('target'); cd('target');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',2,'study',0); 
    SPREP_func(subj_num,EEG,p,1,0,'target') % (EEG,p,purpose=sLORETA,alphaOn=0,condition='Target')

    % open Dataset #3 (Probe) -----------------------------------------------
    mkdir('probe'); cd('probe');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',3,'study',0); 
    SPREP_func(subj_num,EEG,p,1,0,'probe') % (EEG,p,purpose=sLORETA,alphaOn=0,condition='Probe')

    % open Dataset #4 (Irrelevant) -----------------------------------------------
    mkdir('irrelevant'); cd('irrelevant');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',4,'study',0); 
    SPREP_func(subj_num,EEG,p,1,0,'irrelevant') % (EEG,p,purpose=sLORETA,alphaOn=0,condition='Irrelevant')
    
    clear savename;
    cd('..'); cd('..'); % go to subj_Num folder
    
    % log
    fprintf(logFile, ['%.2f - Save ERSP sLORETA txt file: On \n'], toc);
else
    % log when not save sLORETA .txt file
    fprintf(logFile, ['%.2f - Save ERSP sLORETA txt file: Off \n'], toc);
end

%% ERSP
% making folders 
mkdir('TF_analysis'); cd('TF_analysis');
mkdir('Normal'); cd('Normal'); mkdir('ERSP data'); mkdir('ITC data');
cd('..');
mkdir('Bootstrap'); cd('Bootstrap'); mkdir('ERSP data'); mkdir('ITC data');
cd('..');

% log
fprintf(logFile, ['%.2f - ERSP process : Time(%s) Cycle(%.1f %.1f) Freq(%.1f %.1f) Baseline(%s)\n'],...
    toc,num2str(p.wavelet.time),p.wavelet.cycle(1),p.wavelet.cycle(2),p.wavelet.freqs(1),p.wavelet.freqs(2),num2str(p.wavelet.baseline));

% open Dataset #2 (Target) -----------------------------------------------
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',2,'study',0); 
eeglab redraw;
% (EEG,p,purpose=ERSP,alphaOn=0,condition='Target')
cd('Normal'); SPREP_func(subj_num,EEG,p,2,0,'Target'); cd('..');
cd('Bootstrap'); SPREP_func(subj_num,EEG,p,2,1,'Target'); cd('..'); 
% save current EEGLAB file as .set file
if p.wavelet.saveeach == 1
cwd = pwd; name_save = [subj_num,'_CIT_ERSP_tar.set'];
EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
% log
fprintf(logFile, ['%.2f - Save ERSP .set file: %s \n'], toc, name_save);
end

% open Dataset #3 (Probe) -----------------------------------------------
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'retrieve',3,'study',0); 
eeglab redraw;
% (EEG,p,purpose=ERSP,alphaOn=0,condition='Probe')
cd('Normal'); SPREP_func(subj_num,EEG,p,2,0,'Probe'); cd('..');
cd('Bootstrap'); SPREP_func(subj_num,EEG,p,2,1,'Probe'); cd('..'); 
% save current EEGLAB file as .set file
if p.wavelet.saveeach == 1
cwd = pwd; name_save = [subj_num,'_CIT_ERSP_pro.set'];
EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
% log
fprintf(logFile, ['%.2f - Save ERSP .set file: %s \n'], toc, name_save);
end

% open Dataset #4 (Irrelevant)--------------------------------------------
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 3,'retrieve',4,'study',0); 
eeglab redraw;
% (EEG,p,purpose=ERSP,alphaOn=0,condition='Irrelevant')
cd('Normal'); SPREP_func(subj_num,EEG,p,2,0,'Irrelevant'); cd('..');
cd('Bootstrap'); SPREP_func(subj_num,EEG,p,2,1,'Irrelevant'); cd('..'); 
% save current EEGLAB file as .set file
if p.wavelet.saveeach == 1
cwd = pwd; name_save = [subj_num,'_CIT_ERSP_irr.set'];
EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
eeglab redraw;
% log
fprintf(logFile, ['%.2f - Save ERSP .set file: %s \n'], toc, name_save);
end

fprintf(logFile, ['\n']);
cd('..');
% variables which stored are matrix consist of frequency X times
% upperside is lower frequencies, downside is higher frequencies.
% following 'Frequency' matrix.
% ITC matrix which stored are also same.

end % if p.wavelet.calcOn == 1

%% Re-epoch current datasets for ERP
% It means Epoch range to -200 ~ 1500 for ERP
ALLEEG = pop_delset( ALLEEG, [2  3  4] ); eeglab redraw;

% low-pass filter to cut-off frequency at 30Hz
EEG = pop_eegfiltnew(EEG, [], p.ERP.lowpassfilter); % edge 27 for 30 cut-off
EEG = eeg_checkset( EEG ); eeglab redraw;
fprintf(logFile, ['%.2f - ERP Low-passed %dHz (30 cut-off)  \n'], toc,p.ERP.lowpassfilter);

% Re-epoch to short range (-200~1500ms)
EEG = pop_select( EEG,'time',p.ERP.time );
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 
eeglab redraw;
fprintf(logFile, ['%.2f - ERP re-epoch to short size: %.1f ~ %.1f \n'], toc,p.ERP.time(1),p.ERP.time(2));

% Baseline correction to zero-mean (-200ms~0ms)
EEG = pop_rmbase( EEG, p.ERP.baseline);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 
eeglab redraw;
fprintf(logFile, ['%.2f - ERP baseline removal : %d ~ %d \n'], toc,p.ERP.baseline(1),p.ERP.baseline(2));

% Copy current dataset 3 to 4
[ALLEEG EEG CURRENTSET] = pop_copyset( ALLEEG, 3, 4);
eeglab redraw;

EEG = pop_selectevent( EEG, 'type',Index_tr,'renametype','Target','deleteevents','off','deleteepochs','off','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_pr,'renametype','Probe','deleteevents','off','deleteepochs','off','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_ir ,'renametype','Irrelevant','deleteevents','off','deleteepochs','off','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'overwrite','on','gui','off'); 
EEG = eeg_checkset( EEG );
eeglab redraw;

% Save
mkdir('ERP'); cd('ERP');
cwd = pwd; name_save = [subj_num,'_CIT_ERP.set'];
EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
cd('..');
fprintf(logFile, ['%.2f - Save ERP .set file: %s \n'], toc, name_save);

ALLEEG = pop_delset( ALLEEG, [4] ); % remove 4th file(re-naming whole ERP dataset)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'retrieve',3,'study',0); % go 3rd file(rmremove, no re-naming)
eeglab redraw;

% Re-divide files by triggers (target/probe/irrelevant)
EEG = pop_selectevent( EEG, 'type',Index_tr,'deleteevents','on','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Target','gui','off'); Trial_tr = EEG.trials;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',3,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_pr,'deleteevents','on','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Probe','gui','off'); Trial_pr = EEG.trials;
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'retrieve',3,'study',0); 
EEG = eeg_checkset( EEG );
EEG = pop_selectevent( EEG, 'type',Index_ir,'deleteevents','on','deleteepochs','on','invertepochs','off');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname','Irrelevant','gui','off'); Trial_ir = EEG.trials;
eeglab redraw;

%log
fprintf(logFile, ['%.2f - Divide ERP file to each condition: Target(%d TR) / Probe(%d TR) / Irrelevant(%d TR) \n'], toc, Trial_tr, Trial_pr, Trial_ir);
clear Trial_tr Trial_pr Trial_ir

%% save ERP epoch sLORETA files
if p.sLORETA == 1
    cd('sLORETA');
    mkdir('ERP'); cd('ERP');
    % open Dataset #2 (Target) -----------------------------------------------
    mkdir('target'); cd('target');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'retrieve',4,'study',0); 
    SPREP_func(subj_num,EEG,p,1,0,'target') % (EEG,p,purpose=sLORETA,alphaOn=0,condition='Target')

    % open Dataset #3 (Probe) -----------------------------------------------
    mkdir('probe'); cd('probe');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',5,'study',0); 
    SPREP_func(subj_num,EEG,p,1,0,'probe') % (EEG,p,purpose=sLORETA,alphaOn=0,condition='Probe')

    % open Dataset #4 (Irrelevant) -----------------------------------------------
    mkdir('irrelevant'); cd('irrelevant');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'retrieve',6,'study',0); 
    SPREP_func(subj_num,EEG,p,1,0,'irrelevant') % (EEG,p,purpose=sLORETA,alphaOn=0,condition='Irrelevant')

    clear savename;
    cd('..'); cd('..'); % go to subj_Num folder
    fprintf(logFile, ['%.2f - Save ERP sLORETA txt file: On \n'], toc);
else
    fprintf(logFile, ['%.2f - Save ERP sLORETA txt file: Off \n'], toc);
end

%% plotting and save ERP
mkdir('ERP'); cd('ERP');
mkdir('ERP channels'); cd('ERP channels');

% plot comparison ERP
for i = 1:EEG.nbchan-4
pop_comperp( ALLEEG, 1, [4:6] ,[],'addavg','off','addstd','off','addall','on','subavg','on','diffavg','on','diffstd','off','chans',i,'alpha',0.05,'lowpass',30,'tplotopt',{'title' EEG.chanlocs(i).labels 'ydir' 1});
set(gcf,'Position',[50 50 800 600]);
legend('Target(b)','Probe(r)','Irrelevant(g)','location','southwest');
savename = sprintf('[%s]ERP_%s.png',datestr(now,'yyyymmdd'),EEG.chanlocs(i).labels);
%export_fig(savename,gcf); close figure 2;
pic = getframe(gcf);
imwrite(pic.cdata,savename); close figure 2; 
end
cd('..'); %go to 'subj_num \ ERP' folder

% open Dataset #2 (Target)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'retrieve',4,'study',0); 
eeglab redraw;
% plot scalp ERP (Target)
if p.ERP.multiERP == 1
    figure; pop_plottopo(EEG, 1:EEG.nbchan-4, 'Target', 0, 'ydir',1,'ylim',[-15 15] );
    set(gcf,'Position',[50 50 1024 768]);
    savename = sprintf('[%s]Scalp ERP_target.png',datestr(now,'yyyymmdd'));
    pic = getframe(gcf);
    imwrite(pic.cdata,savename); close figure 2;
end
% save current EEGLAB file as .set file
if p.ERP.saveeach == 1
    cwd = pwd; name_save = [subj_num,'_CIT_ERP_tar.set'];
    EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    % log
    fprintf(logFile, ['%.2f - Save ERP .set file: %s \n'], toc, name_save);
end

% open Dataset #3 (Probe)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',5,'study',0); 
eeglab redraw;
% plot scalp ERP (Probe)
if p.ERP.multiERP == 1
    figure; pop_plottopo(EEG, 1:EEG.nbchan-4, 'Probe', 0, 'ydir',1,'ylim',[-15 15] );
    set(gcf,'Position',[50 50 1024 768]);
    savename = sprintf('[%s]Scalp ERP_Probe.png',datestr(now,'yyyymmdd'));
    pic = getframe(gcf);
    imwrite(pic.cdata,savename); close figure 2;
end
% save current EEGLAB file as .set file
if p.ERP.saveeach == 1
    cwd = pwd; name_save = [subj_num,'_CIT_ERP_pro.set'];
    EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    % log
    fprintf(logFile, ['%.2f - Save ERP .set file: %s \n'], toc, name_save);
end

% open Dataset #4 (Irrelevant)
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'retrieve',6,'study',0); 
eeglab redraw;
% plot scalp ERP (Irrelevant)
if p.ERP.multiERP == 1
    figure; pop_plottopo(EEG, 1:EEG.nbchan-4, 'Irrelevant', 0, 'ydir',1,'ylim',[-15 15] );
    set(gcf,'Position',[50 50 1024 768]);
    savename = sprintf('[%s]Scalp ERP_Irrelevant.png',datestr(now,'yyyymmdd'));
    pic = getframe(gcf);
    imwrite(pic.cdata,savename); close figure 2;
end
% save current EEGLAB file as .set file
if p.ERP.saveeach == 1
    cwd = pwd; name_save = [subj_num,'_CIT_ERP_irr.set'];
    EEG = pop_saveset( EEG, 'filename',name_save,'filepath',cwd);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    eeglab redraw;
    % log
    fprintf(logFile, ['%.2f - Save ERP .set file: %s \n'], toc, name_save);
end


% ERP data save to .mat file -------------------------------------------
if p.ERP.savemat == 1
    mkdir('ERP data'); cd('ERP data');
    % save target ERP data
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'retrieve',4,'study',0); eeglab redraw;
    v_ERP = mean(EEG.data,3);
    savename = sprintf('[%s]ERP data_target',datestr(now,'yyyymmdd'));
    save(savename, 'v_ERP'); v_ERP=[];
    % save probe ERP data
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',5,'study',0); eeglab redraw;
    v_ERP = mean(EEG.data,3);
    savename = sprintf('[%s]ERP data_probe',datestr(now,'yyyymmdd'));
    save(savename, 'v_ERP'); v_ERP=[];
    % save irrelevant ERP data
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'retrieve',6,'study',0); eeglab redraw;
    v_ERP = mean(EEG.data,3);
    savename = sprintf('[%s]ERP data_irrelevant',datestr(now,'yyyymmdd'));
    save(savename, 'v_ERP'); v_ERP=[];
    cd('..')
    % log
    fprintf(logFile, ['%.2f - Save ERP .mat files: On \n'], toc);
else
    fprintf(logFile, ['%.2f - Save ERP .mat files: Off \n'], toc);
end

%% Remove monopolar EOG channels, and plot topograph
if p.ERP.topo == 1
    % making folder
    mkdir('ERP topo'); cd('ERP topo');
    % Target topo
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'retrieve',4,'study',0); eeglab redraw;
    SPREP_func(subj_num,EEG,p,3,[],'target'); % function 3 is ERP topofunction
    % Probe topo
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'retrieve',5,'study',0); eeglab redraw; 
    SPREP_func(subj_num,EEG,p,3,[],'probe');
    % Irrelevant topo
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'retrieve',6,'study',0); eeglab redraw; 
    SPREP_func(subj_num,EEG,p,3,[],'irrelevant');
    cd('..');
    % log
    fprintf(logFile, ['%.2f - Save ERP topo files: On \n'], toc);
else
    fprintf(logFile, ['%.2f - Save ERP topo files: Off \n'], toc);
end

cd('..'); cd('..'); % go to 'EEG' folder
%% end message
fprintf(logFile, ['\n \n SPREP _ ERSP & ERP is done \n \n']);
fclose(logFile);

disp('======================================== ');
disp('  SPREP_ERSP&ERP Preprocessing is done');
disp('======================================== ');

close all; clear EEG ALLCOM ALLEEG ALLERP ALLERPCOM CURRENTERP CURRENTSET;
clear CURRENTSTUDY eeglabUpdater ERP LASTCOM plotset PLUGINLIST STUDY;
clear subj_num subj_group subj_target subj_probe logFile v_ERP savename;
clear Index_tr Index_pr Index_ir pic;
pause(2);
% end for Listnum
end

