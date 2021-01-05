function event = ft_read_event_hcnlab_oldtasks(cfg,data,task)

wrkdir = pwd;

datafldr = fileparts(cfg.dataset);
cd(datafldr);
cd .. % go back 1 step to Neural\ folder to see all Blocks

subject = cfg.subject;
block = cfg.blocks;

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

% count all Block folders
blockDirs = dir('Block*');

% for each Block folder, find any bhv2 files
bhvfiles.name = [];
for b=1:length(blockDirs)
    foo = dir(fullfile(blockDirs(b).name,'*.bhv2'));
    if length(foo)>1
        ft_error('More that one bhv2 file found in Block folder.')
    elseif ~isempty(foo)
        bhvfiles(b).name = foo.name;
    end
end
bhvfiles(arrayfun(@(x) isempty(x.name),bhvfiles)) = [];

% sort all bhv2 files across blocks by datetime
bvhdatetimes = NaT();
for b=1:length(bhvfiles)
    nameinfo = strsplit(bhvfiles(b).name,'_');
    nameinfo = [nameinfo{1} nameinfo{2}];
    t = datetime(nameinfo,'inputformat','yyMMddHHmm');
    bvhdatetimes(b) = t;
end
[~,idx] = sort(bvhdatetimes);
bhvfiles = bhvfiles(idx);
blockDirs = blockDirs(idx);

% create uber behavioral_code vector with all ML codes across bhv files
bhv_data = [];
for b=1:length(bhvfiles)
    eval(['data' num2str(b) ' =mlread(fullfile(blockDirs(b).name,bhvfiles(b).name));']);
    bhv_data = [bhv_data eval(['data' num2str(b)])];
end

bhv_cell = {};
for t=1:length(bhv_data)
    bhv_cell = cat(1,bhv_cell, bhv_data(t).BehavioralCodes.CodeNumbers);
end
bhv_codes = cell2mat(bhv_cell)';
trl_bhv_idx = zeros(length(bhv_cell),2);
s = 1;
bhv_rt = zeros(size(bhv_data,2),1);
for c=1:length(bhv_cell) 
    bhv_rt(c) = bhv_data(c).ReactionTime/1000;
    trl_bhv_idx(c,:) = [c s];
    s = s + length(bhv_cell{c});
end

% find where data from AO matches uber behavioral_code vector
new_ao_codes = nan(1,length(data(2,:)));
if length(data(2,:))>=length(bhv_codes)
    chk = strfind(data(2,:),bhv_codes); % using all of bvh_codes
    trl_st = 1;
else
    chk = strfind(bhv_codes,data(2,:)); % using only part of bhv_codes
    if isempty(chk)
        data_starts = find(data(2,:)==9);
        data_ends = find(data(2,:)==18);
        d_st_en = [];e=1;
        for s=1:length(data_starts)
            if data_ends(e)>data_starts(s) && (s==length(data_starts) || data_ends(e)>data_starts(s+1))
                d_st_en = [d_st_en;data_starts(s) data_ends(e)];
                e=e+1;
            else
                if s==length(data_starts)
                    d_st_en = [d_st_en;data_starts(s) size(data,2)];
                else
                    d_st_en = [d_st_en;data_starts(s) data_starts(s+1)-1];
                end
            end
        end
        trl_ns = cell(size(d_st_en,1),2);
        for d=1:size(d_st_en,1)
            trlcodes = data(2,d_st_en(d,1):d_st_en(d,2));
            trl_ns{d,1} = trl_bhv_idx(ismember(trl_bhv_idx(:,2),strfind(bhv_codes,trlcodes)','rows'),1);
            trl_ns{d,2} = trlcodes';
            if any(ismember(trlcodes,[6 8]))
                stim_idx = find(trlcodes==4);
                resp_idx = find(trlcodes==6 | trlcodes==8);
                rt_difs = abs((data(1,d_st_en(d,1)+resp_idx-1) - data(1,d_st_en(d,1)+stim_idx-1))/44000 - bhv_rt);
                trl_ns{d,1} = find(rt_difs<0.01);
            end
        end
        done = all(cellfun(@length,trl_ns(:,1))==1);
        while ~done
            for d=1:size(d_st_en,1)
                if length(trl_ns{d,1})>1
                    if length(trl_ns{d,1})==size(trl_bhv_idx,1)
                        if d-1>0
                            trl_ns{d,1} = trl_ns{d-1,1}+1;
                        end
                    else
                        if d-1>0
                            bef_trn = trl_ns{d-1,1};
                        else
                            if length(trl_ns{d+1,1})==1
                                bef_trn = trl_ns{d+1,1}-2;
                            else
                                bef_trn = -inf;
                            end
                        end
                        if d+1<=size(d_st_en,1)
                            aft_trn = trl_ns{d+1,1};
                        else
                            if length(trl_ns{d-1,1})==1
                                aft_trn = trl_ns{d-1,1}+2;
                            else
                                aft_trn = inf;
                            end
                        end
                        if length(bef_trn)==1 && length(aft_trn)==1
                            cur_pos_trls = trl_ns{d,1}>bef_trn;
                            cur_pos_trls = cur_pos_trls & ...
                                trl_ns{d,1}<aft_trn;
                            cur_pos_trls = trl_ns{d,1}(cur_pos_trls);
                            if length(cur_pos_trls)>1
                                if any(ismember(cur_pos_trls,bef_trn+1)) && ...
                                        all(~ismember(cur_pos_trls,aft_trn-1))
                                    cur_pos_trls = bef_trn+1;
                                elseif any(ismember(cur_pos_trls,aft_trn-1)) && ...
                                        all(~ismember(cur_pos_trls,bef_trn+1))
                                    cur_pos_trls = aft_trn-1;
                                end
                            end
                        elseif length(bef_trn)==1
                            cur_pos_trls = trl_ns{d,1}(trl_ns{d,1}>bef_trn);
                            if any(ismember(cur_pos_trls,bef_trn+1)) && any(ismember(aft_trn,bef_trn+2))
                                cur_pos_trls = bef_trn+1;
                                if d+1<=size(d_st_en,1)
                                    trl_ns{d+1,1} = bef_trn+2;
                                end
                            end
                        elseif length(aft_trn)==1
                            cur_pos_trls = trl_ns{d,1}(trl_ns{d,1}<aft_trn);
                            if any(ismember(cur_pos_trls,aft_trn-1)) && any(ismember(bef_trn,aft_trn-2))
                                cur_pos_trls = aft_trn-1;
                                if d-1>0
                                    trl_ns{d-1,1} = aft_trn-2;
                                end
                            end
                        else
                            cur_pos_trls = trl_ns{d,1};
                        end
                        trl_ns{d,1} = cur_pos_trls;
                    end % if length(trl_ns{d,1})==size(trl_bhv_idx,1)
                end % if length(trl_ns{d,1})>1
            end % for d=1:size(d_st_en,1)
            done = all(cellfun(@length,trl_ns(:,1))==1);
        end % while ~done
        trls = cell2mat(trl_ns(:,1));
        trldif = diff(trls);
        blkchng = [1;find(trldif~=1)];
        conc_blks = nan(size(trls,1),max(cellfun(@length,trl_ns(:,2))));
        for b=2:length(blkchng)
            conc_blks(blkchng(b-1):blkchng(b),1) = b-1;
        end
        conc_blks(blkchng(b)+1:end,1) = b;
        for t=1:size(trl_ns,1)
            conc_blks(t,1:length(trl_ns{t,2})) = repmat(conc_blks(t,1),[1 length(trl_ns{t,2})]);
        end
        conc_blks = conc_blks';
        conc_blks = conc_blks(:);
        conc_blks(isnan(conc_blks)) = [];
        bhv_data = bhv_data(trls);
        bhv_codes = cell2mat(trl_ns(:,2))';
        chk=1;
    else
        trl_st = find(trl_bhv_idx(:,2)==chk);
    end % if isempty(chk)
end % if length(data(2,:))>=length(bhv_codes)

if isempty(chk)
    ft_warning(['Datetimes of bhv files do not match order of behavioral codes '...
        'in AO data file. Will attempt to find matching codes using different ' ...
        'permutations of the order of the bhv files.']);
    
    opts = perms(1:length(blockDirs));
    opts(ismember(opts,idx,'rows'),:) = [];
    c=1;
    while isempty(chk) && c<=size(opts,1)
        dels = repmat('\b',1,23);
        fprintf([dels 'Checking permutation %d'],c);
        
        bhv_data = [];
        for b=opts(c,:)
            eval(['data' num2str(b) ' =mlread(fullfile(blockDirs(b).name,bhvfiles(b).name));']);
            bhv_data = [bhv_data eval(['data' num2str(b)])];
        end
        bhv_cell = {};
        for t=1:length(bhv_data)
            bhv_cell = cat(1,bhv_cell, bhv_data(t).BehavioralCodes.CodeNumbers);
        end
        bhv_codes = cell2mat(bhv_cell)';
        trl_bhv_idx = zeros(length(bhv_cell),2);
        s = 1;
        for c=1:length(bhv_cell)
            trl_bhv_idx(c,:) = [c s];
            s = s + length(bhv_cell{c});
        end
        
        if length(data(2,:))>=length(bhv_codes)
            chk = strfind(data(2,:),bhv_codes);
            trl_st = 1;
        else
            chk = strfind(bhv_codes,data(2,:));
            trl_st = find(trl_bhv_idx(:,2)==chk);
        end
        c=c+1;
    end
    if c>size(opts,1)
        ft_error(['No permutation of behavioral codes from ML match behavioral ' ...
            'codes from AO. Please check data.'])
    end
    ft_warning(fprintf(['Found correct permutation of behavioral codes from ML file. '...
        'The correct ordering is:\n %s\n'],bhvfiles(opts(c-1,:)).name));
end

if length(data(2,:))>=length(bhv_codes)
    new_ao_codes(chk:chk+length(bhv_codes)-1) = bhv_codes;
else
    new_ao_codes = data(2,:);
    bhv_codes = bhv_codes(chk:chk+length(new_ao_codes)-1);
    bhv_data = bhv_data(trl_st:end);
end

% do eventtype sorting based on info in bhv files
if any(isnan(new_ao_codes))
    ft_warning(['Warning: data from AO contains more behavioral codes than '...
        'data from ML. Those extra codes from AO will be provided event.type = ' ...
        '''nan''.']);
end

idx = regexp(cfg.dataset,'Block(s)?\d\w*/');
block = str2double(cfg.dataset(idx+5:find(cfg.dataset=='/',1,'last')-1));
if isnan(block)
    if exist('conc_blks','var')
        block = conc_blks';
    else
        block = ones(size(new_ao_codes));
    end
else
    block = repmat(block,size(new_ao_codes));
end

event.type=[]; event.sample=[]; event.value=[]; event.reactiontime = []; event.task = []; event.block = []; event.trial = []; event.depth = []; event.region = [];
prvtp = 'nan';
tr = 0;
switch task
    case {'SimonTask_v8';'SimonTask_v9'}
        for a=1:length(new_ao_codes)
            tp = 'nan';
            if new_ao_codes(a)==9
                tr=tr+1;
                
                curCond = cell2mat(fields(bhv_data(tr).TaskObject.CurrentConditionInfo));
                codes = bhv_data(tr).BehavioralCodes.CodeNumbers';
                rt = bhv_data(tr).ReactionTime;
                
                colorchk = contains(curCond,'blue','ignorecase',1);
                sidechk = contains(curCond,'_R');
                if colorchk && sidechk % congruent, should press right
                    truth = 'Cong_Right';
                elseif colorchk && ~sidechk % incongruent, should press right
                    truth = 'Incong_Right';
                elseif ~colorchk && sidechk % incongruent, should press left
                    truth = 'Incong_Left';
                else
                    truth = 'Cong_Left';
                end
                
                if ~any(ismember([6 8],codes))
                    resp = '_NoResp';
                    rtval = [];
                elseif ismember(6,codes)
                    resp = '_Left';
                    rtval = 6;
                else
                    resp = '_Right';
                    rtval = 8;
                end
                tp = [truth resp];
                if find(codes==6)<find(codes==4)
                    tp = 'Abort_Early_Left';
                elseif find(codes==8)<find(codes==4)
                    tp = 'Abort_Early_Right';
                end
            else
                tp = prvtp;
            end
            event(a).type = tp;
            event(a).value = data(2,a);
            event(a).sample = data(1,a);
            event(a).trial = tr;
            event(a).block = block(a);
            event(a).depth = depth;
            event(a).region = region;
            event(a).reactiontime = rt;
            event(a).rtvalue = rtval;
            prvtp = tp;
        end
    case 'SimonTask_v10'
        for a=1:length(new_ao_codes)
            tp = 'nan';
            if new_ao_codes(a)==9
                tr=tr+1;
                
                curCond = bhv_data(tr).TaskObject.CurrentConditionInfo.Stim2;
                codes = bhv_data(tr).BehavioralCodes.CodeNumbers';
                rt = bhv_data(tr).ReactionTime;
                
                colorchk = contains(curCond,'blue','ignorecase',1);
                sidechk = contains(curCond,'_R');
                if colorchk && sidechk % congruent, should press right
                    truth = 'Cong_Right';
                elseif colorchk && ~sidechk % incongruent, should press right
                    truth = 'Incong_Right';
                elseif ~colorchk && sidechk % incongruent, should press left
                    truth = 'Incong_Left';
                else
                    truth = 'Cong_Left';
                end
                
                if ~any(ismember([6 8],codes))
                    resp = '_NoResp';
                    rtval = [];
                elseif ismember(6,codes)
                    resp = '_Left';
                    rtval = 6;
                else
                    resp = '_Right';
                    rtval = 8;
                end
                tp = [truth resp];
                if find(codes==6)<find(codes==4)
                    tp = 'Abort_Early_Left';
                elseif find(codes==8)<find(codes==4)
                    tp = 'Abort_Early_Right';
                end
            else
                tp = prvtp;
            end
            event(a).type = tp;
            event(a).value = data(2,a);
            event(a).sample = data(1,a);
            event(a).trial = tr;
            event(a).block = block(a);
            event(a).depth = depth;
            event(a).region = region;
            event(a).reactiontime = rt;
            event(a).rtvalue = rtval;
            prvtp = tp;
        end
    case {'StopTask_v2';'StopTask_v3'}
        for a=1:length(new_ao_codes)
            tp = 'nan';
            if new_ao_codes(a)==9
                tr=tr+1;
                
                curCond = bhv_data(tr).TaskObject.CurrentConditionInfo.Stim3;
                codes = bhv_data(tr).BehavioralCodes.CodeNumbers';
                rt = bhv_data(tr).ReactionTime;
                
                gochk = contains(curCond,'go','ignorecase',1);
                sidechk = contains(curCond,'_left','ignorecase',1);
                if gochk && sidechk % go trial, left arrow
                    truth = 'Go_Left';
                elseif gochk && ~sidechk % go trial, right arrow
                    truth = 'Go_Right';
                elseif ~gochk && sidechk % stop trial, left arrow
                    truth = 'Stop_Left';
                else
                    truth = 'Stop_Right';
                end
                
                if find(codes==16)<find(codes==5)
                    resp = '_LeftBeforeStop';
                    rtval = 16;
                elseif find(codes==17)<find(codes==5)
                    resp = '_RightBeforeStop';
                    rtval = 17;
                elseif~any(ismember([10 12],codes))
                    resp = '_NoResp';
                elseif ismember(10,codes)
                    resp = '_Left';
                    rtval = 10;
                elseif ismember(12,codes)
                    resp = '_Right';
                    rtval = 12;
                elseif ismember(11,codes)
                    resp = '_Left_ERROR';
                    rtval = 11;
                elseif ismember(13,codes)
                    resp = '_Right_ERROR';
                    rtval = 13;
                end
                tp = [truth resp];
                if find(codes==14)<find(codes==3)
                    tp = 'Abort_Early_Left';
                    rtval = 14;
                elseif find(codes==15)<find(codes==3)
                    tp = 'Abort_Early_Right';
                    rtval = 15;
                end
            else
                tp = prvtp;
            end
            event(a).type = tp;
            event(a).value = data(2,a);
            event(a).sample = data(1,a);
            event(a).trial = tr;
            event(a).block = block(a);
            event(a).depth = depth;
            event(a).region = region;
            event(a).reactiontime = rt;
            event(a).rtvalue = rtval;
            prvtp = tp;
        end
end
cd(wrkdir);

end