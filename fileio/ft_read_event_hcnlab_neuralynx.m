function [event] = ft_read_event_hcnlab_neuralynx(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

hdr          = ft_getopt(varargin, 'header');
cfg          = ft_getopt(varargin,'cfg');
if isempty(cfg)
    headerformat = ft_getopt(varargin, 'headerformat');
    eventformat  = ft_getopt(varargin, 'eventformat');
    dataformat   = ft_getopt(varargin, 'dataformat');
    task   = ft_getopt(varargin, 'task');
else
    headerformat = cfg.headerformat;
    eventformat = cfg.eventformat;
    dataformat = cfg.dataformat;
    task = cfg.trialdef.task;
end
if ~contains(task,'task','ignorecase',true)
    task = [task 'Task'];
else
    task(strfind(task,'ask')-1) = upper(task(strfind(task,'ask')-1));
end
task(1) = upper(task(1)); 

[task_event_dict,task] = get_event_types(task);
old_tasks = {'SimonTask_v8';'SimonTask_v9';'SimonTask_v10';'StopTask_v2';'StopTask_v3'};

[data] = ft_read_data(cfg.dataset,'header',hdr);
if contains(cfg.dataset,'ncs')
    thresh1 = -3000;
    pls_idx = data<thresh1;
    pls1_st = find(diff([0 pls_idx])==1)-2;
    thresh2 = 3000;
    pls_idx = data>thresh2;
    pls2_st = find(diff([0 pls_idx])==1)-2;
    
    pulses = [pls1_st pls2_st;ones(size(pls1_st)) 2*ones(size(pls2_st))];
    [~,idx] = sort(pulses(1,:));
    pulses = pulses(:,idx);
    data = pulses;
else
    ft_error('ft_read_event_hcnlab_neuralynx.m not completed for reading non-ncs events'); 
end

event.type=[]; event.sample=[]; event.value=[]; event.task = []; event.block = []; event.trial = []; event.depth = []; event.region = [];

idx = regexp(cfg.dataset,'Block[0-9]+/');
block = str2double(cfg.dataset(idx+5:find(cfg.dataset=='/',1,'last')-1));

type = [];
if ismember(task,old_tasks)
    event = ft_read_event_hcnlab_oldtasks(cfg,data,task);
else
    parts = strsplit(cfg.dataset,'/');
    subject = parts{5};
    file = parts{end};
    ft_warning('ft_read_event_hcnlab_neuralynx.m not complete to identify recording depth and location information.');
    %%% NEED A DIFFERENT METHOD TO IDENTIFY DEPTH  & REGION INFORMATION
%     depthsstr = '[dD][\w\W]+[fF]';
%     [s,e] = regexp(file,depthsstr);
%     depth = str2double(file(s+1:e-1));
%     
%     recLocMetaDataFile = '/media/HCN_Data/SubjectMetaData/CSV_files/RecLocsBySubjectTask/metaData.csv';
%     recLoc_tbl = readtable(recLocMetaDataFile);
%     idx_s = strcmpi(recLoc_tbl.Subject,subject);
%     idx_d = recLoc_tbl.MicroDepth==depth;
%     
%     idx = idx_s & idx_d;
%     if ~any(idx)
%         ft_error(['Subject ' subject ' not found in RecLocsBySubjectTask/metaData.csv. Please add.']);
%     elseif sum(idx)>1
%         ft_error(['Combination of subject ' subject ' and depth ' depth ' occurs more than once '...
%             'in RecLocsBySubjectTask/metaData.csv. Please fix.']);
%     end
%     region = recLoc_tbl.Region(idx);
    %%%
    depth = []; region = [];
    
    if contains(cfg.dataset,'ncs')
        t=0;
        for j=1:size(data,2)
            event(j).value = data(2,j);
            event(j).sample = data(1,j);
            type_val = data(2,j);
            idx = [task_event_dict{1,:}]==type_val;
            type = task_event_dict{2,idx};
            event(j).type = type;
            event(j).trial = j;
            event(j).block = block;
            event(j).depth = depth;
            event(j).region = region;
        end
    else
        ft_error('ft_read_event_hcnlab_neuralynx.m not completed for reading non-ncs events');
%         t=0;
%         for j=1:size(data,2)
%             event(j).value = data(2,j);
%             event(j).sample = data(1,j);
%             if data(2,j)==9
%                 t=t+1;
%                 type_val = data(2,j+find(ismember(data(2,j+1:end),[task_event_dict{1,:}]),1,'first'));
%                 idx = [task_event_dict{1,:}]==type_val;
%                 type = task_event_dict{2,idx};
%             end
%             event(j).type = type;
%             event(j).trial = t;
%             event(j).block = block;
%             event(j).depth = depth;
%             event(j).region = region;
%         end
    end
end

[event(:).task] = deal(task);

if ~isempty(event)
  % make sure that all required elements are present
  if ~isfield(event, 'type'),     ft_error('type field not defined for each event');     end
  if ~isfield(event, 'sample'),   ft_error('sample field not defined for each event');   end
  if ~isfield(event, 'value'),    for i=1:length(event), event(i).value = [];    end; end
  if ~isfield(event, 'offset'),   for i=1:length(event), event(i).offset = [];   end; end
  if ~isfield(event, 'duration'), for i=1:length(event), event(i).duration = []; end; end
end

% make sure that all numeric values are double
if ~isempty(event)
  for i=1:length(event)
    if isnumeric(event(i).value)
      event(i).value = double(event(i).value);
    end
    event(i).sample    = double(event(i).sample);
    event(i).offset    = double(event(i).offset);
    event(i).duration  = double(event(i).duration);
  end
end

if ~isempty(event)
  % sort the events on the sample on which they occur
  % this has the side effect that events without a sample number are discarded
  sample = [event.sample];
  if ~all(isnan(sample))
    [~, indx] = sort(sample);
    event = event(indx);
  end
end

% apply the optional filters
event = ft_filter_event(event, varargin{:}); %% need to check this to see if it
                                              % will work for my needs

if isempty(event)
  % ensure that it has the correct fields, even if it is empty
  event = struct('type', {}, 'value', {}, 'sample', {}, 'offset', {}, 'duration', {});
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION that maps event codes to event types for each task
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [event_type_dict,task] = get_event_types(task)
expsDir = fullfile('/media/HCN_Data/Experiments');
func1 = @(x) contains(x.name,task);
func2 = @(y) str2double(y.name(strfind(y.name,'_v')+2:end));

verchk = arrayfun(@(x) real(str2double(x)),task);
if ~any(verchk)
    ft_warning(['No version number identified for task ' task '. Looking in '...
        'task directory to identify latest version.']);
%     taskDir = dir(fullfile(expsDir,task));
%     verDirs = taskDir(arrayfun(func1,taskDir));
%     vers = arrayfun(func2,verDirs);
%     latestVer = vers(vers==max(vers));
%     task = [task '_v' num2str(latestVer)];
    task = getLatestTaskVersion(task);
else
    if ~contains(task,'_v')
        idx = find(isnan(verchk),1,'last');
        task = [task(1:idx) '_v' task(idx+1:end)];
    end
end

switch task
    case {'SimonTask_v8';'SimonTask_v9';'SimonTask_v10'}
        event_type_dict = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,...
            24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45;...
            'Start_Trial','Fix_On','Fix_Off','Cue_On','Cue_Off','LEFT_Button_On','LEFT_Button_Off',...
            'RIGHT_Button_On','ML_Start_Trial','RIGHT_Button_Off','Early_Left_On_pre_recording',...
            'Early_Left_Off_pre_recording','Early_Right_On_pre_recording','Early_Right_Off_pre_recording',...
            'Early_Left_On_fixation','Early_Left_Off_fixation','Early_Right_On_fixation','ML_End_Trial',...
            'Early_Right_Off_fixation','Temp20','EarlyPress_ERROR_pre_recording','EarlyPress_ERROR_fixation',...
            'EarlyPress_ERROR_Pic_On','EarlyPress_ERROR_Pic_Off','NoResponse_ERROR','NoResponse_ERROR_Pic_On',...
            'NoResponse_ERROR_Pic_Off','Temp28','Temp29','Temp30','LateResponse_Left_On',...
            'LateResponse_Left_Off','LateResponse_Right_On','LateResponse_Right_Off',...
            'LateResponse_ERROR','Temp36','Temp37','Temp38','Temp39','Temp40','pre_recording_time_on',...
            'pre_recording_time_off','post_recording_time_on','post_recording_time_off','End_Trial'};
    case {'SimonTask_v11';'SimonTask_v12'} % need to complete latest version that ports all necessary trial codes to AO
        event_type_dict = {10,11,12,13,14,15,16,17,19,20,21,22,23,24;...
            'Cong_Left_Left','Cong_Left_Right','Cong_Left_NoResp','Cong_Right_Left',...
            'Cong_Right_Right','Cong_Right_NoResp','Abort_Early_Left','Abort_Early_Right',...
            'Incong_Left_Left','Incong_Left_Right','Incong_Left_NoResp','Incong_Right_Left',...
            'Incong_Right_Right','Incong_Right_NoResp'};
    case {'StopTask_v2';'StopTask_v3'}
        event_type_dict = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21;...
            'Start_Trial','Fix_On','Fix_Off','Cue1_On','Cue1_Off','Cue2_On','Cue2_Off',...
            'NoEvent1','ML_Start_Trial','Left_Button_On','Left_Button_On_ERROR',...
            'Right_Button_On','Right_Button_On_ERROR','Early_Left_Button_On_fixation',...
            'Early_Right_Button_On_fixation','Early_Left_Button_On_cue1','Early_Right_Button_On_cue1',...
            'ML_End_Trial','NoResponse','NoResponse_ERROR','End_Trial'};
    case {'StopTask_v4'} % need to complete latest version that ports all necessary trial codes to AO
        event_type_dict = {12,13,14,15,16,17,19,20,21,22,23,24,25,26,27,28,29,30;...
            'Go_Left_Left','Go_Left_Right','Go_Left_NoResp','Go_Right_Left',...
            'Go_Right_Right','Go_Right_NoResp','Stop_Left_Left','Stop_Left_Right',...
            'Stop_Left_NoResp','Stop_Right_Left','Stop_Right_Right','Stop_Right_NoResp',...
            'Abort_Early_Left','Abort_Early_Right','Stop_Left_LeftBeforeStop',...
            'Stop_Left_RightBeforeStop','Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
    case 'ChoiceTask_v1'
        event_type_dict = {1;...
            'Cue_On'};
    case 'BDMAuctionTask_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuralynx.'])
    case 'RatingTask_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuralynx.'])
    case 'RL4Task_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuralynx.'])
    case 'FaceWordLocalizerTask_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuralynx.'])
    otherwise
        ft_error(['Task ' task ' not recognized. Please check or add to ' ...
            'ft_read_event_hcnlab_neuralynx.'])
end
end