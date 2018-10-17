function varargout = process_conn_PSI( varargin )
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
    sProcess.Comment     = 'Phase Slope Index (PSI)';
    sProcess.FileTag     = 'PSI';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Connectivity_Fieldtrip';
    sProcess.Index       = 996; % add arbitrary number. it is not important.
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data',     'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
    sProcess.isSeparator = 1; % draw seperate line in GUI menus
    
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
    
    % === FREQ DECOMPOSITION METHOD
    sProcess.options.label3.Comment = '<BR><U><B>Estimator options</B></U>:';
    sProcess.options.label3.Type    = 'label';
    sProcess.options.tapers.Comment = {'Hanning window', 'Multitaper window'};
    sProcess.options.tapers.Type    = 'radio';
    sProcess.options.tapers.Value   = 2;
    %sProcess.options.tfcyclemin.Comment = 'Set Min Cycles (Wavelet):';
    %sProcess.options.tfcyclemin.Type    = 'value';
    %sProcess.options.tfcyclemin.Value   = {3,'cycles',2};
    %sProcess.options.tfcyclemax.Comment = 'Set Max Cycles (Wavelet):';
    %sProcess.options.tfcyclemax.Type    = 'value';
    %sProcess.options.tfcyclemax.Value   = {10,'cycles',2};
    sProcess.options.freqsmoothf.Comment = 'Freqeuncy Smoothing Factor(only Multitaper):';
    sProcess.options.freqsmoothf.Type    = 'value';
    sProcess.options.freqsmoothf.Value   = {3,'¡¾Hz',2};
    
    % === FREQUENCY RESOLUTION
    sProcess.options.label4.Comment = '';
    sProcess.options.label4.Type    = 'label';
    %sProcess.options.freqres.Comment = 'Frequency resolution (Interval):';
    %sProcess.options.freqres.Type    = 'value';
    %sProcess.options.freqres.Value   = {2,'Hz',2};
    % === HIGHEST FREQUENCY OF INTEREST
    sProcess.options.maxfreq.Comment = 'Highest frequency of interest:';
    sProcess.options.maxfreq.Type    = 'value';
    sProcess.options.maxfreq.Value   = {100,'Hz',2};
    % === LOWEST FREQUENCY OF INTEREST
    sProcess.options.lowfreq.Comment = 'Lowest frequency of interest:';
    sProcess.options.lowfreq.Type    = 'value';
    sProcess.options.lowfreq.Value   = {2,'Hz',2};
    
    % === SELECT BANDWIDTH
    sProcess.options.label5.Comment = '<BR><U><B>Band Width</B></U>:';
    sProcess.options.label5.Type    = 'label';
    sProcess.options.usingbandwidth.Comment = {'Rayleigh frequency(Fieldtrip default, 1/T)','Input Value'};
    sProcess.options.usingbandwidth.Type    = 'radio';
    sProcess.options.usingbandwidth.Value   = 1;
    sProcess.options.bandwidth.Comment = 'Bandwidth size (Input Value selected only):';
    sProcess.options.bandwidth.Type    = 'value';
    sProcess.options.bandwidth.Value   = {2,'¡¾Hz',2};
    
    % === OUTPUT MODE    
    sProcess.options.label6.Comment = '<BR><U><B>Output configuration</B></U>:';
    sProcess.options.label6.Type    = 'label';
    sProcess.options.outputmode.Comment = {'Fieldtrip recommend connectivity analysis on multiple trials (one file)'};
    sProcess.options.outputmode.Type    = 'radio';
    sProcess.options.outputmode.Value   = 1;

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
    
    OPTIONS.Method = 'wpli';
    OPTIONS.RemoveEvoked = sProcess.options.removeevoked.Value;
    
    
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
    
    clear Num_data 
    
    cfg            = [];
    Data           = ft_preprocessing(cfg,Data);
    
    % Re-epoch to time of interest
    cfg            = [];
    cfg.toilim     = [sProcess.options.timewindow.Value{1}(1) sProcess.options.timewindow.Value{1}(2)];
    Data = ft_redefinetrial(cfg, Data);
    
    % Frequency decomposition
    cfg            = [];
    cfg.output     = 'fourier'; % return complex Fourier-spectra
    cfg.keeptrials = 'yes';
    %cfg.foi        = [sProcess.options.lowfreq.Value{1} : sProcess.options.freqres.Value{1} : sProcess.options.maxfreq.Value{1}];
    %cfg.toi        = OPTIONS.TimeWindow(1):(1000/Data.fsample)/1000:OPTIONS.TimeWindow(2);
    cfg.foilim     = [sProcess.options.lowfreq.Value{1} sProcess.options.maxfreq.Value{1}];
    
    if sProcess.options.removeevoked.Value == 1
        cfg.demean = 'yes';
    else
        cfg.demean = 'no';
    end
    
    if sProcess.options.tapers.Value == 1
        cfg.taper  = 'hanning';
    elseif sProcess.options.tapers.Value == 2
        cfg.taper  = 'dpss';
        cfg.tapsmofrq = sProcess.options.freqsmoothf.Value{1};
    end
        
        
        cfg.method = 'mtmfft';
        
        %cfg.width  = linspace(sProcess.options.tfcyclemin.Value{1},sProcess.options.tfcyclemax.Value{1},length(cfg.foi));
    
    freqdata       = ft_freqanalysis(cfg, Data);
    

    % Calculation PSI connectivity
    cfg            = [];
    cfg.method     = 'psi';
    if sProcess.options.usingbandwidth.Value == 1
        cfg.bandwidth = 1./(OPTIONS.TimeWindow(2) - OPTIONS.TimeWindow(1));
    elseif sProcess.options.usingbandwidth.Value == 2
        cfg.bandwidth = sProcess.options.bandwidth.Value{1}.*2;
    end
    psimatrix        = ft_connectivityanalysis(cfg, freqdata);

    % extract result data and change format to brainstorm file
    TMP_result = psimatrix.psispctrm;

    %TMP_result = permute(TMP_result,[1,2,4,3]); %permute change dimension order..!
    result = reshape(TMP_result, Num_scout^2, 1, size(TMP_result,3));
    
    ResultMat = db_template('timefreqmat');
    ResultMat.TF = result; % save calculated matrix file
    ResultMat.Comment      = sprintf(['PSI(%.0fms-%.0fms)'],sProcess.options.timewindow.Value{1}(1)*1000, sProcess.options.timewindow.Value{1}(2)*1000);
    ResultMat.DataType     = 'matrix';
    ResultMat.DataFile     = [];
    ResultMat.Time         = [sProcess.options.timewindow.Value{1}(1), sProcess.options.timewindow.Value{1}(2)];
    ResultMat.TimeBands{1} = 'psi'; ResultMat.TimeBands{2} = ResultMat.Time(1); ResultMat.TimeBands{3} = ResultMat.Time(2);
    ResultMat.Freqs        = psimatrix.freq;
    ResultMat.RefRowNames  = psimatrix.label;
    ResultMat.RowNames     = psimatrix.label;
    ResultMat.Measure      = 'other';
    ResultMat.Method       = 'psi'; 
    ResultMat.SurfaceFile  = temp.SurfaceFile;
    ResultMat.Atlas        = temp.Atlas;
    ResultMat.nAvg         = length(sInputs);
    ResultMat.Options      = OPTIONS;
    %ResultMat.Options.Freqs = [sProcess.options.lowfreq.Value{1} : sProcess.options.freqres.Value{1} : sProcess.options.maxfreq.Value{1}];
    ResultMat.Options.Freqs = psimatrix.freq;
    ResultMat.Options.isSymmetric = logical(0);
    ResultMat.Options.Outputmode   = 'concat';
    ResultMat.Options.isScoutA     = logical(0);
    ResultMat.Options.isScoutB     = logical(0);
    ResultMat.Options.sScoutA      = [];
    ResultMat.Options.sScoutB      = [];
    ResultMat.Options.iOutputStudy = sInputs(1).iStudy;

    clear TMP_result Data freqdata psimatrix 
    

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
