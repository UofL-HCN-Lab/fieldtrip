function event = ft_get_events_from_monkeylogic(ao_data,cfg,task_event_dict,fs)

if ~exist('fs','var') || isempty(fs)
    hdr = ft_read_header(cfg.headerfile, 'headerformat', cfg.headerformat,'chantype',cfg.chantype);
    fs = hdr.portFs;
end
wrkdir = pwd;

datafldr = fileparts(cfg.dataset);
cd(datafldr);
cd .. % go back 1 step to Neural\ folder to see all Blocks

subject = cfg.subject;
block = cfg.blocks;

% get metadata to add to events struct later on
recLocMetaDataFile = fullfile(cfg.metaDataDir,'metaData.csv');
opts = detectImportOptions(recLocMetaDataFile);
opts = setvaropts(opts,{'Subject','Block'},'type','string');
recLoc_tbl = readtable(recLocMetaDataFile,opts);
idx_s = strcmpi(recLoc_tbl.Subject,subject);
blocks = cellfun(@str2double,cellfun(@strsplit,recLoc_tbl.Block,'uniformoutput',0),'uniformoutput',0);
idx_b = cellfun(@(x) isequal(block,x),blocks);
idx_t = strcmpi(recLoc_tbl.Task,cfg.task);
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
elseif contains(cfg.chantype,'macro')
    depth = recLoc_tbl.MacroDepth(idx);
elseif contains(cfg.chantype,'eeg')
    depth = recLoc_tbl.SEEG{idx};
end

if ~any(idx)
    ft_error(['Subject ' subject ' not found in RecLocsBySubjectTask/metaData.csv. Please add.']);
elseif sum(idx)>1
    ft_error(['Combination of subject ' subject ' and depth ' num2str(depth') ' occurs more than once '...
        'in RecLocsBySubjectTask/metaData.csv. Please fix.']);
end
region = recLoc_tbl.Region(idx);

% load bhv file
datafldr = fileparts(cfg.dataset);
cd(datafldr);

bhvfiles = dir('*.bhv2');
if length(bhvfiles)>1
    ft_error('More that one bhv2 file found in Block folder.')
end

bhvdata = mlread(bhvfiles(1).name);

% get all bhv codes
bhvcodes = [];
bhvtimes = [];
for t=1:length(bhvdata)
    bhvcodes = [bhvcodes;bhvdata(t).BehavioralCodes.CodeNumbers];
    bhvtimes = [bhvtimes;bhvdata(t).AbsoluteTrialStartTime + bhvdata(t).BehavioralCodes.CodeTimes];
end
bhvcodes = bhvcodes';
bhvtimes = bhvtimes';

% mapping from 'incorrect' eventcode in AO to correct associated code in
% bhv file

if isequal(cfg.eventvalue,[1 16])
    bhv_eventvalue = [9 18]; % for trial start
else
    ft_error(['Event value ' num2str(cfg.eventvalue) ' not added to mapping in ft_get_events_from_monkeylogic.m']);   
end

% compare # stim on events AO to ML to check they are consistent
% also, get relative times between repeats of eventvalue codes (for
% comparison)

% first, ML
if length(bhv_eventvalue)>1
    idxs = false(length(bhv_eventvalue),length(bhvcodes)+length(bhv_eventvalue)-1);
    start = length(bhv_eventvalue)+1;
    for e=1:length(bhv_eventvalue)
        idxs(e,start-e:(length(bhvcodes)+start-e-1)) = bhvcodes==bhv_eventvalue(e);
    end
    idx_ml = all(idxs,1); idx_ml = idx_ml(1:(end-length(length(bhv_eventvalue))));
    mltimes = [bhvtimes(1,find(idx_ml)-1);bhvtimes(1,find(idx_ml))];
    mltimes = diff(mltimes,1,1);
    
    foo = num2str(bhvcodes);
    nevents_ml = length(strfind(foo,num2str(bhv_eventvalue)));
else
    idx_ml = bhvcodes==bhv_eventvalue;
    mltimes = diff(bhvtimes(idx_ml));
    
    nevents_ml = sum(idx_ml);
end

% next, AO
if length(cfg.eventvalue)>1
    idxs = false(length(cfg.eventvalue),size(ao_data,2)+length(cfg.eventvalue)-1);
    start = length(cfg.eventvalue)+1;
    for e=1:length(cfg.eventvalue)
        idxs(e,start-e:(size(ao_data,2)+start-e-1)) = ao_data(2,:)==cfg.eventvalue(e);
    end
    idx_ao = all(idxs,1); idx_ao = idx_ao(1:(end-length(length(cfg.eventvalue))));
    aotimes = [ao_data(1,find(idx_ao)-1);ao_data(1,find(idx_ao))];
    aotimes = diff(aotimes,1,1)*1000/fs;
    
    foo = num2str(ao_data(2,:));
    nevents_ao = length(strfind(foo,num2str(cfg.eventvalue)));
else
    idx_ao = ao_data(2,:)==cfg.eventvalue;
    aotimes = diff(ao_data(1,idx_ao))*1000/fs;
    
    nevents_ao = sum(idx_ao);
end
idx_ao = find(idx_ao);

% uh oh, number of events don't match. attempt parsimonious comparison by
% relative durations between successive events
if nevents_ao~=nevents_ml
    ft_warning(['Number of events from ML bhv2 files does not match number of events from AO. '...
        'Attempting to determine best mapping by comparing relative code times. Please check '...
        'summary stats before accepting.']);
    
    if nevents_ao<nevents_ml
        ao_shorter = true; longer = mltimes; shorter = aotimes;
        ttl_txt = 'ML shifted by ';
    else
        ao_shorter = false; longer = aotimes; shorter = mltimes;
        ttl_txt = 'AO shifted by ';
    end
    
    % compare durration differences by shifting data with longer length one
    % event at a time to find the shift that produces a comparison that is
    % optimal (most duration differences < 1 ms)
    tries = length(longer)-length(shorter)+1;
    comps = zeros(tries,length(shorter));
    for t=1:tries
        comps(t,:) = longer(t:(length(shorter)+t-1)) - shorter;
    end
    best = find(sum(abs(comps)<=1,2)==max(sum(abs(comps)<=1,2)));
    
    % if more than one option is identified as "best"
    if length(best)>1
        f=figure;
        for b=1:length(best)
            subplot(1,length(best),b);
            histogram(comps(best(b),:));
            title({['Best Auto Option ' num2str(b)];[ttl_txt num2str(best(b)-1)]});
        end
        prompt = ['Using the plotted histograms, please identify which option number\n'...
            'gives the most parsimonious alignment between ML and AO trial data.\n'...
            'Input nan or empty vector ([]) if no alignment matches.\n'];
        resp = input(prompt);
        close(f);
        
        % no options are good, error out
        if isempty(resp) || isnan(resp)
            ft_error('No alignment acceptable between ML and AO trial data. Please verify correct eventvalues for ML and AO.');
        end
        best = best(resp);
    else
        % otherwise, confirm best, or manual choose
        f=figure;
        histogram(comps(best,:));
        title({'Best Auto Option';[ttl_txt num2str(best-1)]});
        
        % ask for user confirmation
        prompt = ['Auto-detected best option displayed in histogram. Accept (y/n)?\n'];
        resp = input(prompt,'s');
        close(f);
        % user rejected auto-detection, allow for manual comparison
        if ~strncmpi('yes',resp,length(resp))
            fprintf(['Auto-detected best option rejected. Plotting all possible options for manual\n' ...
                'selection.\n']);
            
            f=figure;
            for t=1:tries
                subplot(1,tries,t);
                histogram(comps(t,:));
                title({['Option ' num2str(t)];[ttl_txt num2str(t-1)]});
            end
            prompt = ['Using the plotted histograms, please identify which option number\n'...
                'gives the most parsimonious alignment between ML and AO trial data.\n'...
                'Input nan or empty vector ([]) if no alignment matches.\n'];
            best = input(prompt);
            close(f);
            
            % no options are good, error out
            if isempty(best) || isnan(best)
                ft_error('No alignment acceptable between ML and AO trial data. Please verify correct eventvalues for ML and AO.');
            end
        end
    end
    
    % if ao has fewer events, remove extra trials from bhvdata, otherwise
    % remove "best minus 1"-most events found in AO data
    if ao_shorter
        bhvdata = bhvdata(best:(length(shorter)+best-1));
    else
        idx_ao = idx_ao(best:(length(shorter)+best-1));
    end
end
aosamples = ao_data(1,idx_ao);

event.task = []; event.block = []; event.trial = []; event.type=[]; event.sample=[]; event.value=[]; event.reactiontime = [];  event.rtvalue = []; event.depth = []; event.region = [];
switch lower([cfg.task '_v' cfg.taskver])
    case {'simontask_v13','simontaskstim_v13'}
        j=1;
        t=0;
        for tr=1:length(bhvdata)
            trcodes = bhvdata(tr).BehavioralCodes.CodeNumbers';
            trtimes = bhvdata(tr).BehavioralCodes.CodeTimes';
            
            % find idx0 (alignment between matched AO and ML event codes),
            % and set trtimes to be relative to trtimes(idx0)
            idx0 = false(length(bhv_eventvalue),length(trcodes)+length(bhv_eventvalue)-1);
            start = length(bhv_eventvalue)+1;
            for e=1:length(bhv_eventvalue)
                idx0(e,start-e:(length(trcodes)+start-e-1)) = trcodes==bhv_eventvalue(e);
            end
            idx0 = all(idx0,1); idx0 = idx0(1:(end-length(length(bhv_eventvalue))));
            trtimes = trtimes-trtimes(idx0);
            
            % convert trtimes to relative sample number from matched AO
            % event sample
            trsamples = fix(trtimes/1000*fs);
            for e=find(idx0):length(trcodes)
                if e==find(idx0)
                    event(j).value = trcodes(1);
                else
                    event(j).value = trcodes(e);
                end

                if event(j).value==9
                    stim_idx = e+find(ismember(trcodes(e+1:end),4),1,'first');
                    abort_idx = e+find(ismember(trcodes(e+1:end),[15 16 17]),1,'first');
                    if isempty(stim_idx) & ~isempty(abort_idx) % trial aborted before stimulus, skip
                        skip=true;
                    else
                        skip=false;
                        t=t+1;
                    end
                    
                    type_val = trcodes(e+find(ismember(trcodes(e+1:end),[task_event_dict{1,:}]),1,'first'));
                    idx = [task_event_dict{1,:}]==type_val;
                    type = task_event_dict{2,idx};
                    
                    if isempty(stim_idx) & isempty(abort_idx) & any(contains(task_event_dict(2,:),'Instructions'))
                        stim_idx = e+find(ismember(trcodes(e+1:end),8),1,'first');
                    end
                    butt_idx = e+find(ismember(trcodes(e+1:end),[6 7]),1,'first');
                    rt = round((trsamples(butt_idx)-trsamples(stim_idx))*1000/fs); % for rt in ms
                    if isempty(rt), rt = nan; end
                    rtval = trcodes(butt_idx);
                end
                if ~skip
                    event(j).sample = trsamples(e) + aosamples(tr);
                    event(j).trial = t;
                    event(j).type = type;
                    event(j).block = str2double(regexprep(num2str(block),' ',''));
                    event(j).depth = depth;
                    event(j).region = region;
                    event(j).reactiontime = rt;
                    event(j).rtvalue = rtval;
                    j=j+1;
                end
            end
        end
        
    otherwise
        ft_error(['ft_get_events_from_monkeylogic.m not coded to work with '...
            [cfg.task '_v' cfg.taskver] '. Please check task is correct in cfg '...
            'or update code.']);
end

cd(wrkdir);