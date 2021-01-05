function [event] = ft_read_event_hcnlab_neuroomega(varargin)
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

[task_event_dict,task] = get_event_types_local(task);
old_tasks = {'SimonTask_v8';'SimonTask_v9';'SimonTask_v10';'StopTask_v2';'StopTask_v3'};

fields_orig=who(hdr.orig); %getting digital event channels
trig_chan = fields_orig{strcmp(fields_orig,'CPORT__1')};

event.task = []; event.block = []; event.trial = []; event.type=[]; event.sample=[]; event.value=[]; event.reactiontime = [];  event.rtvalue = []; event.depth = []; event.region = [];
data = hdr.orig.(trig_chan);

% remove partial trials at beginning of data (specifically for 20190313_JH
% because recording was started in the middle of a previous trial)
frst_fullTrial_idx = find(data(2,:)==9,1,'first');
if ~isempty(frst_fullTrial_idx)
    data = data(:,frst_fullTrial_idx:end);
    good_AO_codes = true;
else
    % if no 9 exists, this must be from subject data when an incorrect
    % behavioral code mapping was inputted into AO due to a splitter box
    % from ML not working as intended. we will assume the data starts from
    % the beginning (e.g., trial 1)
    good_AO_codes = false;
end

block = cfg.blocks;
% [idxb, idxe] = regexp(cfg.dataset,'Block(s)?[0-9]+\w*/');
% blockdata = strsplit(cfg.dataset(idxb:(idxe-1)),'s');
% if length(blockdata)==1
%     block = str2double(blockdata{1}((strfind(blockdata{1},'k')+1):end));
% else
%     blockdata = blockdata{2};
%     blockdata = strsplit(blockdata,'and');
%     block = str2double(blockdata);
% end

type = [];
if ismember(task,old_tasks)
    data(:,data(2,:)==27) = []; % removes spurious 27 event code
    event = ft_read_event_hcnlab_oldtasks(cfg,data,task);
else
    subject = cfg.subject;
    
    recLocMetaDataFile = fullfile(cfg.metaDataDir,'metaData.csv');
    opts = detectImportOptions(recLocMetaDataFile);
    opts = setvaropts(opts,{'Subject','Block'},'type','string');
    recLoc_tbl = readtable(recLocMetaDataFile,opts);
    idx_s = strcmpi(recLoc_tbl.Subject,subject);
    idx_t = strcmpi(recLoc_tbl.Task,cfg.task);
    blocks = cellfun(@str2double,cellfun(@strsplit,recLoc_tbl.Block,'uniformoutput',0),'uniformoutput',0);
    idx_b = cellfun(@(x) isequal(block,x),blocks);
    idx = idx_s & idx_b & idx_t;
    hems = recLoc_tbl.Hemisphere(idx);
    if length(unique(hems))>1
        if any(contains(cfg.channel,'_01___'))
            idx_h = strcmp(recLoc_tbl.Hemisphere,'Left');
        else
            idx_h = strcmp(recLoc_tbl.Hemisphere,'Right');
        end
        idx = idx & idx_h;
    end
    if contains(cfg.chantype,'micro')
        depth = recLoc_tbl.MicroDepth(idx);
        idx_d = ~isnan(recLoc_tbl.MicroDepth);
    elseif contains(cfg.chantype,'macro')
        depth = recLoc_tbl.MacroDepth(idx);
        idx_d = ~isnan(recLoc_tbl.MacroDepth);
    elseif contains(cfg.chantype,'eeg')
        depth = recLoc_tbl.SEEG(idx);
        func = @(x) str2double(strsplit(x));
        depth = cellfun(func,depth,'uniformoutput',false);
        idx_d = ~cellfun(@isempty,recLoc_tbl.SEEG);
    end
    idx = idx & idx_d;
    if iscell(depth)
        func = @(x) any(isnan(x));
        depth(cellfun(func,depth)) = [];
    else
        depth(isnan(depth)) = [];
    end

    if ~any(idx)
        ft_error(['Subject ' subject ' not found in RecLocsBySubjectTask/metaData.csv. Please add.']);
    elseif sum(idx)>1
        ft_error(['Combination of subject ' subject ' and depth ' num2str(depth(:)') ' occurs more than once '...
            'in RecLocsBySubjectTask/metaData.csv. Please fix.']);
    end
    region = recLoc_tbl.Region(idx);
    
    if ~good_AO_codes
        event = ft_get_events_from_monkeylogic(data,cfg,task_event_dict);
    else
        t=0;
        skip=false;
        for j=1:size(data,2)
            event(j).sample = nan;
            if data(2,j)==9
                trl_end_idx = find(data(2,j:end)==18,1,'first');
                stim_idx = j+find(ismember(data(2,j+1:j+trl_end_idx-1),4),1,'first');
                abort_idx = j+find(ismember(data(2,j+1:j+trl_end_idx-1),[15 16 17]),1,'first');
                if isempty(stim_idx) & ~isempty(abort_idx) % trial aborted before stimulus, skip
                    skip=true;
                else
                    skip=false;
                    t=t+1;
                end
                
                type_val = data(2,j+find(ismember(data(2,j+1:j+trl_end_idx-1),[task_event_dict{1,:}]),1,'first'));
                idx = [task_event_dict{1,:}]==type_val;
                type = task_event_dict{2,idx};
                
                
                if isempty(stim_idx) & isempty(abort_idx) & any(contains(task_event_dict(2,:),'Instructions'))
                    stim_idx = j+find(ismember(data(2,j+1:j+trl_end_idx-1),8),1,'first');
                end
                butt_idx = j+find(ismember(data(2,j+1:j+trl_end_idx-1),[6 7]),1,'first');
                rt = round((data(1,butt_idx)-data(1,stim_idx))*1000/hdr.portFs); % for rt in ms
                if isempty(rt), rt = nan; end
                rtval = data(2,butt_idx);
            end
            if ~skip
                event(j).value = data(2,j);
                event(j).sample = data(1,j);
                event(j).type = type;
                event(j).trial = t;
                event(j).block = str2double(regexprep(num2str(block),' ',''));
                event(j).depth = depth;
                event(j).region = region;
                event(j).reactiontime = rt;
                event(j).rtvalue = rtval;
            end
        end
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
function [event_type_dict,task] = get_event_types_local(task)

verchk = arrayfun(@(x) real(str2double(x)),task);
if ~any(verchk)
    ft_warning(['No version number identified for task ' task '. Looking in '...
        'task directory to identify latest version.']);
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
    case {'SimonTask_v11';'SimonTask_v12'}
        event_type_dict = {10,11,12,13,14,15,16,17,19,20,21,22,23,24;...
            'Cong_Left_Left','Cong_Left_Right','Cong_Left_NoResp','Cong_Right_Left',...
            'Cong_Right_Right','Cong_Right_NoResp','Abort_Early_Left','Abort_Early_Right',...
            'Incong_Left_Left','Incong_Left_Right','Incong_Left_NoResp','Incong_Right_Left',...
            'Incong_Right_Right','Incong_Right_NoResp'};
    case {'SimonTask_v13'}
        event_type_dict = {8,10,11,12,13,14,15,16,17,19,20,21,22,23,24;...
            'Instructions','Cong_Left_Left','Cong_Left_Right','Cong_Left_NoResp','Cong_Right_Left',...
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
    case 'BDMAuctionTask_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega.'])
    case 'RatingTask_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega.'])
    case 'RL4Task_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega.'])
    case 'FaceWordLocalizerTask_v1'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega.'])
    otherwise
        ft_error(['Task ' task ' not recognized. Please check or add to ' ...
            'ft_read_event_hcnlab_neuroomega.'])
end
end