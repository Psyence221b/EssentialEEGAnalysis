% Change trigger code for 2020 Inhibition-Induced Forgetting experiment
% HPLAB. Feb, 2020. written by SIU

%% Parameter
p.Numblock = 6;      % number of blocks
p.Numtrial = 62;     % number of trials within each block
p.Nummaxstim = 5;    % maximum number of stimulus in each trial

%% Remove first 'boundary' event if exist
if strcmp(EEG.event(1).type,  'boundary') == 1, EEG = pop_editeventvals(EEG,'delete',1); end
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

%% Load and re-allocation of event structure
% from EEG structure to matrix form
if isnumeric(EEG.event(1).type) == 1 && isnumeric(EEG.event(2).type) == 1
    for i = 1:size(EEG.event,2), Event_matrix(1,i) = EEG.event(i).type; end
elseif ischar(EEG.event(1).type) == 1 && ischar(EEG.event(2).type) == 1
    event_cell = cell(1,size(EEG.event,2));
    for i = 1:size(event_cell,2), event_cell{i} = str2num(EEG.event(i).type); end
    Event_matrix = cell2mat(event_cell);
else
    error('Cannot examine data type of EEG.evet structure');
end

% find specific value by repeated number
Ind_block = find(Event_matrix == 1); % find number 1 as first trial in each block
if size(Ind_block,2) ~= p.Numblock, error('Block number is not matched. Check again!'); end
Ind_block = [Ind_block size(Event_matrix,2)+1];

for i = 1:p.Numblock % block
    Ind_current = Ind_block(i);
    for j = 1:p.Numtrial - 1 % maximum value of trigger in a specific block except the last trial in each block

        % In normal case : 5 event codes in one trial
        if Event_matrix(Ind_current) == j && Event_matrix(Ind_current+3) == 250 && Event_matrix(Ind_current+5) == j+1 
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current:Ind_current+4)]; count_plus = 5;
        % In normal case : 4 event codes in one trial (no response trials)  
        elseif Event_matrix(Ind_current) == j && Event_matrix(Ind_current+3) == 255 && Event_matrix(Ind_current+4) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current:Ind_current+2) 0 Event_matrix(Ind_current+3)]; count_plus = 4;

        % In abnormal case : one entire trial (5 events) missing
        elseif Event_matrix(Ind_current) == j && Event_matrix(Ind_current+3) == 250 && Event_matrix(Ind_current+5) == j+2 
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current:Ind_current+4)]; count_plus = 5;
        % In abnormal case : one entire trial (4 events) missing    
        elseif Event_matrix(Ind_current) == j && Event_matrix(Ind_current+3) == 255 && Event_matrix(Ind_current+4) == j+2 
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current:Ind_current+2) 0 Event_matrix(Ind_current+3)]; count_plus = 4;
        % In abnormal case : code 230 (second) is missing
        elseif Event_matrix(Ind_current) == j && Event_matrix(Ind_current+1) ~= 230
            if Event_matrix(Ind_current+1) == 235 && Event_matrix(Ind_current+4) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current) 0 Event_matrix(Ind_current+1:Ind_current+3)]; count_plus = 4;
            elseif Event_matrix(Ind_current+1) == 235 && Event_matrix(Ind_current+3) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current) 0 Event_matrix(Ind_current+1:Ind_current+2)]; count_plus = 3;
            % case: both code 230 and 235 are missing
            elseif Event_matrix(Ind_current+1) == 250 && Event_matrix(Ind_current+3) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current) 0 0 Event_matrix(Ind_current+1:Ind_current+2)]; count_plus = 3;    
            elseif Event_matrix(Ind_current+1) == 255 && Event_matrix(Ind_current+2) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current) 0 0 0 Event_matrix(Ind_current+1)]; count_plus = 2;    
            end
        % In abnormal case : code 235 (third) is missing
        elseif Event_matrix(Ind_current) == j && Event_matrix(Ind_current+2) ~= 235
            if Event_matrix(Ind_current+2) == 250 && Event_matrix(Ind_current+3) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current:Ind_current+1) 0 Event_matrix(Ind_current+2:Ind_current+3)]; count_plus = 4;    
            elseif Event_matrix(Ind_current+2) == 255 && Event_matrix(Ind_current+2) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = [Event_matrix(Ind_current:Ind_current+1) 0 0 Event_matrix(Ind_current+2)]; count_plus = 3;        
            end
        % In abnormal case : others
        elseif Event_matrix(Ind_current) ~= j && Event_matrix(Ind_current) == j+1
            Mat_reshape((i-1)*p.Numtrial+j,:) = zeros(1,p.Nummaxstim); count_plus = 0;
        end
        Ind_current = Ind_current + count_plus;
    end
    
    % the last trial in each block
    if Event_matrix(Ind_current+3) == 250
    Mat_reshape((i-1)*p.Numtrial+p.Numtrial,:) = Event_matrix(Ind_current:Ind_block(i+1)-1);
    elseif Event_matrix(Ind_current+3) == 255
    Mat_reshape((i-1)*p.Numtrial+p.Numtrial,:) = [Event_matrix(Ind_current:Ind_current+2) 0 Ind_block(i+1)-1];
    end

end

%% Check modified trigger is same with the original data
TMP_check = permute(Mat_reshape, [2,1]);
TMP_check2 = reshape(TMP_check, [1,p.Numblock * p.Numtrial * p.Nummaxstim]);
TMP_check2(TMP_check2 == 0) = [];
TMP_check3 = Event_matrix == TMP_check2;
check = find(TMP_check3 == 0);
if isempty(check) ~= 1, error('Some error occur while change the triggers!'); end

% Check the size of trigger
if size(Mat_reshape,1) ~= p.Numblock * p.Numtrial, error('Some triggers are missing!'); end

%% Change the values of trigger
Mat_change = Mat_reshape;

Ind_zero = find(Mat_change(:,1) == 0); % find where whole trial is missing

Mat_change(:,1) = 1:p.Numblock * p.Numtrial; % change triggers to continuous values across blocks

if isempty(Ind_zero) ~= 1
    for i = 1:length(Ind_zero)
    Mat_change(Ind_zero(i),:) = zeros(1,p.Nummaxstim); % fill with zero(0) the missing trial
    end
end

% find location of values of 0
TMP_Mat_change = reshape(Mat_change',1,[]);
Ind_null = find(TMP_Mat_change == 0);

% change to cell array (for string values)
Mat_change = num2cell(Mat_change);

for i = 1:size(Mat_change,1)
    if Mat_change{i,2} ~= 0, Mat_change{i,2} = [num2str(Mat_change{i,1}) '-sti']; end
    if Mat_change{i,3} ~= 0, Mat_change{i,3} = [num2str(Mat_change{i,1}) '-bep']; end
    if Mat_change{i,4} ~= 0, Mat_change{i,4} = [num2str(Mat_change{i,1}) '-res']; end
    if Mat_change{i,5} ~= 0, Mat_change{i,5} = [num2str(Mat_change{i,1}) '-fed']; end
    if Mat_change{i,1} ~= 0, Mat_change{i,1} = [num2str(Mat_change{i,1}) '-fix']; end
end

% make to one line array
Mat_change = reshape(Mat_change',1,[]);
Mat_change(Ind_null) = [];

if size(Mat_change,2) ~= size(EEG.event,2)
    error('Some error occur while change the triggers!');
end

%% Change EEGLAB structure with the variable 'Mat_change'
for i = 1:size(EEG.event,2)
    EEG.event(i).type = Mat_change{1,i};
end

clear i j count_plus Ind_zero Ind_null Event_cell Event_matrix Ind_block Ind_current 
clear TMP_check TMP_check2 TMP_check3 check TMP_Mat_change

EEG = eeg_checkset( EEG );
eeglab redraw
