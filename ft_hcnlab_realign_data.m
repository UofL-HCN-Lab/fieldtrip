function data = ft_hcnlab_realign_data(data,datatype,eventtype,latency)

if ~exist('latency','var') || isempty(latency)
    latency = [-inf inf];
end
switch datatype
    case {'spike';'spikes';'spk';'spks'}
        switch eventtype
            case {'stim';'stimulus'}
            case {'resp';'response'}
                    for d=1:length(data)
                        all_trials = [data(d).event.trial];
                        if isfield(data(d),'spikeTrials')
                            for c=1:length(data(d).spikeTrials.trial)
                                trls = unique(data(d).spikeTrials.trial{c});
                                for t=1:length(trls)
                                    eidx = all_trials==trls(t);
                                    if ~any(eidx)
                                        disp('hey');
                                    end
                                    rt = data(d).event(eidx).reactiontime/1000; % in s
                                    trlidx = data(d).spikeTrials.trial{c}==trls(t);
                                    data(d).spikeTrials.time{c}(trlidx) = ...
                                        data(d).spikeTrials.time{c}(trlidx)-rt;
                                    
                                    tidx = trlidx & (data(d).spikeTrials.time{c}<latency(1) | ...
                                        data(d).spikeTrials.time{c}>latency(2));
                                    data(d).spikeTrials.time{c}(tidx) = [];
                                    data(d).spikeTrials.trial{c}(tidx) = [];
                                    data(d).spikeTrials.unit{c}(tidx) = [];
                                    data(d).spikeTrials.waveform{c}(:,:,tidx) = [];
                                    data(d).spikeTrials.timestamp{c}(tidx) = [];
                                end
                                data(d).spikeTrials.trialtime = ...
                                    repmat(latency,[size(data(d).spikeTrials.trialtime,1) 1]);
                            end
                        end
                        if isfield(data(d),'spikeTrialsBin')
                            for t=1:length(data(d).spikeTrialsBin.trial)
                                rt = data(d).event(t).reactiontime/1000; % in s
                                data(d).spikeTrialsBin.time{t} = ...
                                    data(d).spikeTrialsBin.time{t}-rt;
                                tidx = data(d).spikeTrialsBin.time{t}>=latency(1) & ...
                                    data(d).spikeTrialsBin.time{t}<=latency(2);
                                data(d).spikeTrialsBin.time{t} = data(d).spikeTrialsBin.time{t}(tidx);
                                data(d).spikeTrialsBin.trial{t} = data(d).spikeTrialsBin.trial{t}(:,tidx);
                            end
                        end
                    end
            otherwise
                ft_error(['Event type ' eventtype ' not yet supported in ft_hcnlab_realign']);
        end
    otherwise
        ft_error(['Data type ' datatype ' not yet supported in ft_hcnlab_realign']);
end