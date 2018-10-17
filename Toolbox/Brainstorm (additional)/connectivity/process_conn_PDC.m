function varargout = process_conn_PDC( varargin )
% PROCESS_ADDCONN: Add connectivity methods by using Fieldtrip functions in
% Brainstorm. add by SIU

% @=============================================================================
% This software is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2013 Brainstorm by the University of Southern California
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Francois Tadel, 2013

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Partial Directed Coherence (PDC)';
    sProcess.FileTag     = 'PDC';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Connectivity_Fieldtrip';
    sProcess.Index       = 1000; % add arbitrary number. it is not important.
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data',     'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
    
    % === CONNECT INPUT
    sProcess = process_corr1n('DefineConnectOptions', sProcess, 1);
    % === REMOVE EVOKED REPONSE
    sProcess.options.removeevoked.Comment = 'Remove evoked response from each trial';
    sProcess.options.removeevoked.Type    = 'checkbox';
    sProcess.options.removeevoked.Value   = 0;
    
    % === SAMPLING RATE
    sProcess.options.label2.Comment = '<BR><U><B>Sampling Frequency</B></U>:';
    sProcess.options.label2.Type    = 'label';
    sProcess.options.sf.Comment = 'Sampling Frequency:';
    sProcess.options.sf.Type    = 'value';
    sProcess.options.sf.Value   = {250,'Hz',2};
    
    % === GRANGER ORDER
    sProcess.options.label3.Comment = '<BR><U><B>Estimator options</B></U>:';
    sProcess.options.label3.Type    = 'label';
    sProcess.options.grangerorder.Comment = 'Multivariate Autoregressive order (default=2):';
    sProcess.options.grangerorder.Type    = 'value';
    sProcess.options.grangerorder.Value   = {2, '', 0};
    % === FREQUENCY RESOLUTION
    sProcess.options.freqres.Comment = 'Frequency resolution (Interval):';
    sProcess.options.freqres.Type    = 'value';
    sProcess.options.freqres.Value   = {2,'Hz',2};
    % === HIGHEST FREQUENCY OF INTEREST
    sProcess.options.maxfreq.Comment = 'Highest frequency of interest:';
    sProcess.options.maxfreq.Type    = 'value';
    sProcess.options.maxfreq.Value   = {100,'Hz',2};
    % === LOWEST FREQUENCY OF INTEREST
    sProcess.options.lowfreq.Comment = 'Lowest frequency of interest:';
    sProcess.options.lowfreq.Type    = 'value';
    sProcess.options.lowfreq.Value   = {2,'Hz',2};
    
    % === OUTPUT MODE    
    sProcess.options.label4.Comment = '<BR><U><B>Output configuration</B></U>:';
    sProcess.options.label4.Type    = 'label';
    sProcess.options.outputmode.Comment = {'Fieldtrip recommend MVAR analysis on multiple trials (one file)'};
    sProcess.options.outputmode.Type    = 'radio';
    sProcess.options.outputmode.Value   = 1;

    sProcess.options.label5.Comment = '<BR><U><B>Result Figure Output configuration</B></U>:';
    sProcess.options.label5.Type    = 'label';    
    % === PICTURE OUTPUT FOLDER
    SelectOptions = {...
        '', ...                            % Filename
        '', ...                            % FileFormat
        'save', ...                        % Dialog type: {open,save}
        'Select output folder...', ...     % Window title
        'ExportData', ...                  % LastUsedDir: {ImportData,ImportChannel,ImportAnat,ExportChannel,ExportData,ExportAnat,ExportProtocol,ExportImage,ExportScript}
        'single', ...                      % Selection mode: {single,multiple}
        'dirs', ...                        % Selection mode: {files,dirs,files_and_dirs}
        {{'.folder'}, 'PNG (*.png)', 'PNG'}}; % Available file formats
    % Option definition
    sProcess.options.outputdir.Comment = 'Figure output folder:';
    sProcess.options.outputdir.Type    = 'filename';
    sProcess.options.outputdir.Value   = SelectOptions;
    % === OUTPUT FILE TAG
    sProcess.options.filetag.Comment = 'Figure output file name:';
    sProcess.options.filetag.Type    = 'text';
    sProcess.options.filetag.Value   = '';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    % Input options
    OPTIONS = process_corr1n('GetConnectOptions', sProcess, sInputs);
    if isempty(OPTIONS)
        OutputFiles = {};
        return
    end
    
    OPTIONS.Method = 'spgranger';
    OPTIONS.RemoveEvoked = sProcess.options.removeevoked.Value;
    OPTIONS.GrangerOrder = sProcess.options.grangerorder.Value{1};
    OPTIONS.MaxFreqRes   = sProcess.options.freqres.Value{1};
    OPTIONS.MaxFreq      = sProcess.options.maxfreq.Value{1};
    
    % check inputs
    if isempty(sProcess.options.outputdir.Value{1}) || isempty(sProcess.options.filetag.Value)
        bst_report('Error', sProcess, [], ['No output folder or name' ]);
        error;
    end
    
    % get data and transform to Fieldtrip format
    Data = [];
    for Num_data = 1:length(sInputs)
        temp = in_bst_data(sInputs(Num_data).FileName);
        Data.trial{Num_data} = temp.Value;
        Data.time{Num_data} = temp.Time;
    end
    Data.fsample = sProcess.options.sf.Value{1};
    Data.label = cell(size(temp.Atlas.Scouts,2),1);
    for Num_scout = 1:size(temp.Atlas.Scouts,2)
    Data.label{Num_scout} = temp.Atlas.Scouts(Num_scout).Label;
    end
    Data.cfg.ntrials = Num_data;
    Data.cfg.triallength = 1;
    Data.cfg.fsample = Data.fsample;
    Data.cfg.nsignal = Num_scout;
    Data.cfg.method = 'ar';
    
    clear Num_data 
    
    % Re-epoch to time of interest
    cfg            = [];
    cfg.toilim     = [sProcess.options.timewindow.Value{1}(1) sProcess.options.timewindow.Value{1}(2)];
    Data = ft_redefinetrial(cfg, Data);
    
    for k = 1:5,     TMP_rank(k) = rank(Data.trial{k}); , end
    fprintf('\nRank of first 5 trials(to check rank dificiency) : %d %d %d %d %d\n',TMP_rank(1),TMP_rank(2),TMP_rank(3),TMP_rank(4),TMP_rank(5));
    fprintf('\nIf error occur in armorf.m, decrease your number of scouts\n');
    
    % Calculate mvar coefficient (multivariate autoregressive coefficient)
    cfg         = [];
    cfg.order   = sProcess.options.grangerorder.Value{1};
    if sProcess.options.removeevoked.Value == 1
        cfg.demean  = 'yes';
    else
        cfg.demean = 'no';
    end
    cfg.toolbox    = 'bsmart';
    mvardata       = ft_mvaranalysis(cfg, Data);
    
    cfg            = [];
    cfg.method     = 'mvar';
    cfg.foi        = [sProcess.options.lowfreq.Value{1} : sProcess.options.freqres.Value{1} : sProcess.options.maxfreq.Value{1}];
    mfreq          = ft_freqanalysis(cfg, mvardata);
    
    % calculation PDC
    cfg            = [];
    cfg.method     = 'pdc';
    granger        = ft_connectivityanalysis(cfg, mfreq);
    
    % reshape for easy to access (by cell)
    PDCresult = cell(Num_scout,Num_scout);
    for row=1:Num_scout
    for col=1:Num_scout
      PDCresult{row,col} = squeeze(granger.pdcspctrm(row,col,:))';
    end
    end    
    
    % plotting & save
    fprintf('\n Preparing for plotting PDC result..... \n');
    cwd = pwd;
    outdir = sProcess.options.outputdir.Value{1};
    outfile = [sProcess.options.filetag.Value '.png'];
    mkdir(outdir);
    cd(outdir)
    
    cfg            = [];
    cfg.parameter  = 'pdcspctrm';
    cfg.zlim       = [0 1];
    ft_connectivityplot(cfg, granger);
    set(gcf,'Position',[10 10 1200 800]);
    pause(2);
        %     figure
        %     for row=1:101
        %     for col=1:10
        %       subplot(10,10,(row-1)*10+col);
        %       plot(granger.freq, squeeze(granger.dtfspctrm(row,col,:)))
        %       ylim([0 1])
        %     end
        %     end
    pic = getframe(gcf);
    imwrite(pic.cdata,outfile); close figure 1;
    
    outfile = ['ft_' sProcess.options.filetag.Value '.mat'];
    save(outfile, 'granger');
    
    outfile = [sProcess.options.filetag.Value '.mat'];
    save(outfile, 'PDCresult');

    cd(cwd); % back to working folder
    
    ResultMat = db_template('timefreqmat');
    ResultMat.Comment      = sprintf(['PDC(%.0fms-%.0fms)'],sProcess.options.timewindow.Value{1}(1)*1000, sProcess.options.timewindow.Value{1}(2)*1000);
    ResultMat.DataType     = 'matrix';
    ResultMat.DataFile     = [];
    ResultMat.Time         = [sProcess.options.timewindow.Value{1}(1), sProcess.options.timewindow.Value{1}(2)];
    ResultMat.TimeBands{1} = 'spgranger'; ResultMat.TimeBands{2} = ResultMat.Time(1); ResultMat.TimeBands{3} = ResultMat.Time(2);
    ResultMat.Freqs        = granger.freq;
    ResultMat.RefRowNames  = granger.label;
    ResultMat.RowNames     = granger.label;
    ResultMat.Measure      = 'other';
    ResultMat.Method       = 'spgranger'; % actually not, but for convinient
    ResultMat.SurfaceFile  = temp.SurfaceFile;
    ResultMat.Atlas        = temp.Atlas;
    ResultMat.nAvg         = length(sInputs);
    ResultMat.Options      = OPTIONS;
    ResultMat.Options.Freqs = [sProcess.options.lowfreq.Value{1} : sProcess.options.freqres.Value{1} : sProcess.options.maxfreq.Value{1}];
    ResultMat.Options.isSymmetric = logical(0);
    ResultMat.Options.Outputmode   = 'concat';
    ResultMat.Options.isScoutA     = logical(0);
    ResultMat.Options.isScoutB     = logical(0);
    ResultMat.Options.sScoutA      = [];
    ResultMat.Options.sScoutB      = [];
    ResultMat.Options.iOutputStudy = sInputs(1).iStudy;
    
    TMP_TF = reshape(PDCresult,1,[]);
    % reshape to Brainstorm spgranger file format
    for i = 1:size(TMP_TF{1},2)
    for j = 1:size(TMP_TF,2)
    ResultMat.TF(j,1,i) = TMP_TF{j}(i);
    end
    end
    
    clear TMP_TF i j row col Data mvardata mfreq granger
    

    % ===== SAVE THE RESULTS =====
    % Get output study
    iOutputStudy = sInputs(1).iStudy;
    sOutputStudy = bst_get('Study', iOutputStudy);

    % Output filename
    fileTag = 'connectn';
    NewFile = bst_process('GetNewFilename', bst_fileparts(sOutputStudy.FileName), ['timefreq_' fileTag '_' OPTIONS.Method]);
    % Save file
    bst_save(NewFile, ResultMat, 'v6');
    % Register in database
    db_add_data(iOutputStudy, NewFile, ResultMat);
    OutputFiles{1} = NewFile;
end





