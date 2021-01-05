function eventtype = ft_trialfun_trial_hcnlab(task,eventtype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION that converts user defined event types into task-defined
% event-type groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tasks = {'SimonTask';'StopTask';'BDMAuctionTask';'RatingTask';'RL4Task';'FaceWordLocalizerTask'};
task = tasks{cellfun(@(x) contains(task,x),tasks)};

switch task
    case {'SimonTask'}
        switch lower(eventtype)
            case 'corr' % all correct
                eventtype = {'Cong_Left_Left','Cong_Right_Right','Incong_Left_Left',...
                    'Incong_Right_Right'};
            case 'incorr' % all incorrect
                eventtype = {'Cong_Left_Right','Cong_Right_Left','Incong_Left_Right',...
                    'Incong_Right_Left'};
            case 'left' % all left stimuli
                eventtype = {'Cong_Left_Left','Cong_Left_Right','Incong_Right_Left',...
                    'Incong_Right_Right'};
            case 'right' % all right stimuli
                eventtype = {'Cong_Right_Left','Cong_Right_Right','Incong_Left_Left',...
                    'Incong_Left_Right'};
            case 'corr_left' % all correct left stimuli
                eventtype = {'Cong_Left_Left','Incong_Right_Right'};
            case 'incorr_left' % all error left stimuli
                eventtype = {'Cong_Left_Right','Incong_Right_Left'};
            case 'corr_right' % all correct right stimuli
                eventtype = {'Cong_Right_Right','Incong_Left_Left'};
            case 'incorr_right' % all error right stimuli
                eventtype = {'Cong_Right_Left','Incong_Left_Right'};
            case 'leftpress' % all left press
                eventtype = {'Cong_Left_Left','Cong_Right_Left','Incong_Left_Left',...
                    'Incong_Right_Left'};
            case 'rightpress' % all right press
                eventtype = {'Cong_Left_Right','Cong_Right_Right','Incong_Left_Right',...
                    'Incong_Right_Right'};
            case 'corr_leftpress' % all correct left press
                eventtype = {'Cong_Left_Left','Incong_Left_Left'};
            case 'incorr_leftpress' % all error left press
                eventtype = {'Incong_Right_Left','Cong_Right_Left'};
            case 'corr_rightpress' % all correct right press
                eventtype = {'Cong_Right_Right','Incong_Right_Right'};
            case 'incorr_rightpress' % all error right press
                eventtype = {'Incong_Left_Right','Cong_Left_Right'};
            case 'corr_cong' % all correct congruent trials
                eventtype = {'Cong_Left_Left','Cong_Right_Right'};
            case 'incorr_cong' % all error congruent trials
                eventtype = {'Cong_Left_Right','Cong_Right_Left'};
            case 'corr_incong' % all correct incongruent trials
                eventtype = {'Incong_Right_Right','Incong_Left_Left'};
            case 'incorr_incong' % all error incongruent trials
                eventtype = {'Incong_Left_Right','Incong_Right_Left'};
            case 'corr_cong_left' % redundant
                eventtype = {'Cong_Left_Left'};
            case 'incorr_cong_left' % redundant
                eventtype = {'Cong_Left_Right'};
            case 'corr_cong_right' % redundant
                eventtype = {'Cong_Right_Right'};
            case 'incorr_cong_right' % redundant
                eventtype = {'Cong_Right_Left'};
            case 'corr_incong_left' % redundant
                eventtype = {'Incong_Right_Right'};
            case 'incorr_incong_left' % redundant
                eventtype = {'Incong_Right_Left'};
            case 'corr_incong_right' % redundant
                eventtype = {'Incong_Left_Left'};
            case 'incorr_incong_right' % redundant
                eventtype = {'Incong_Left_Right'};
            case 'corr_cong_leftpress' % redundant
                eventtype = {'Cong_Left_Left'};
            case 'incorr_cong_leftpress' % redundant
                eventtype = {'Cong_Right_Left'};
            case 'corr_cong_rightpress' % redundant
                eventtype = {'Cong_Right_Right'};
            case 'incorr_cong_rightpress' % redundant
                eventtype = {'Cong_Left_Right'};
            case 'corr_incong_leftpress' % redundant
                eventtype = {'Incong_Left_Left'};
            case 'incorr_incong_leftpress' % redundant
                eventtype = {'Incong_Right_Left'};
            case 'corr_incong_rightpress' % redundant
                eventtype = {'Incong_Right_Right'};
            case 'incorr_incong_rightpress' % redundant
                eventtype = {'Incong_Left_Right'};
            case 'abort' % all early press trials
                eventtype = {'Abort_Early_Left','Abort_Early_Right'};
            case 'noresp' % all no-response trials
                eventtype = {'Cong_Right_NoResp','Cong_Left_NoResp',...
                    'Incong_Right_NoResp','Incong_Left_NoResp'};
            otherwise
                eventtype = {eventtype};
        end
    case {'StopTask'}
        switch lower(eventtype)
            case 'go' % all go trials
                eventtype = {'Go_Left_Left','Go_Left_Right','Go_Left_NoResp',...
                    'Go_Right_Left','Go_Right_Right','Go_Right_NoResp'};
            case 'stop' % all stop trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Left_NoResp',...
                    'Stop_Right_Left','Stop_Right_Right','Stop_Right_NoResp',...
                    'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop',...
                    'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'stop_noearlypress' % exclude stop trials when responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Left_NoResp',...
                    'Stop_Right_Left','Stop_Right_Right','Stop_Right_NoResp'};
            case 'stop_earlypress' % ONLY stop trials when response were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop',...
                    'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'corr' % all correct (go + correct stop trials)
                eventtype = {'Go_Left_Left','Go_Right_Right','Stop_Left_NoResp',...
                    'Stop_Right_NoResp'};
            case 'incorr' % all incorrect (nogo for go trials and go for stop trials)
                eventtype = {'Go_Left_Right','Go_Left_NoResp','Go_Right_Left',...
                    'Go_Right_NoResp','Stop_Left_Left','Stop_Left_Right',...
                    'Stop_Right_Left','Stop_Right_Right','Stop_Left_LeftBeforeStop',...
                    'Stop_Left_RightBeforeStop','Stop_Right_LeftBeforeStop',...
                    'Stop_Right_RightBeforeStop'};
            case 'left' % all left stimulus conditions (correct or incorrect)
                eventtype = {'Go_Left_Left','Go_Left_Right','Go_Left_NoResp',...
                    'Stop_Left_Left','Stop_Left_Right','Stop_Left_NoResp',...
                    'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop','Stop_Left_NoResp'};
            case 'right' % all right stimulus conditions (correct or incorrect)
                eventtype = {'Go_Right_Left','Go_Right_Right','Go_Right_NoResp',...
                    'Stop_Right_Left','Stop_Right_Right','Stop_Right_NoResp',...
                    'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop','Stop_Right_NoResp'};
            case 'leftpress' % all left press (correct or incorrect)
                eventtype = {'Go_Left_Left','Go_Right_Left','Stop_Left_Left',...
                    'Stop_Right_Left','Stop_Left_LeftBeforeStop','Stop_Right_LeftBeforeStop'};
            case 'rightpress' % all right press (correct or incorrect)
                eventtype = {'Go_Left_Right','Go_Right_Right','Stop_Left_Right',...
                    'Stop_Right_Right','Stop_Left_RightBeforeStop','Stop_Right_RightBeforeStop'};
            case 'corr_go' 
                eventtype = {'Go_Left_Left','Go_Right_Right'};
            case 'incorr_go' % all incorrect go trials
                eventtype = {'Go_Left_Right','Go_Right_Left','Go_Left_NoResp','Go_Right_NoResp'};
            case 'incorr_go_resp' % only incorrect go trials with a response
                eventtype = {'Go_Left_Right','Go_Right_Left'};
            case 'incorr_go_noresp' % only incorrect go trials without a response
                eventtype = {'Go_Left_NoResp','Go_Right_NoResp'};
            case 'corr_stop'
                eventtype = {'Stop_Left_NoResp','Stop_Right_NoResp'};
            case 'incorr_stop' % response on stop trials (including responses BEFORE stop signal displays)
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Right_Left',...
                    'Stop_Right_Right','Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop',...
                    'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'incorr_stop_noearlypress' % only responses on stop trials AFTER stop signal displays
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Right_Left',...
                    'Stop_Right_Right'};
            case 'incorr_stop_earlypress' % ONLY responses on stop trials BEFORE stop signal displays
                eventtype = {'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop','Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'corr_leftpress' % all correct left press (go trials only); redundant
                eventtype = {'Go_Left_Left'};
            case 'incorr_leftpress' % all error left press (stop trials only)
                eventtype = {'Go_Right_Left','Stop_Left_Left','Stop_Right_Left',...
                    'Stop_Left_LeftBeforeStop','Stop_Right_LeftBeforeStop'};
            case 'corr_rightpress' % all correct right press (go trials only); redundant
                eventtype = {'Go_Right_Right'};
            case 'incorr_rightpress' % all error right press (stop trials only)
                eventtype = {'Go_Left_Right','Stop_Left_Right','Stop_Right_Right',...
                    'Stop_Left_RightBeforeStop','Stop_Right_RightBeforeStop'};
            case 'go_left' % all go-left trials
                eventtype = {'Go_Left_Left','Go_Left_Right','Go_Left_NoResp'};
            case 'go_right'
                eventtype = {'Go_Right_Left','Go_Right_Right','Go_Right_NoResp'};
            case 'stop_left' % all stop-left trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Left_NoResp',...
                    'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop'};
            case 'stop_left_noearlypress' % exclude stop-left trials when responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Left_NoResp'};
            case 'stop_left_earlypress' % ONLY stop-left trials when response were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop'};
            case 'stop_right' % all stop-right trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Right_Left','Stop_Right_Right','Stop_Right_NoResp',...
                    'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'stop_right_noearlypress' % exclude stop-right trials when responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Right_Left','Stop_Right_Right','Stop_Right_NoResp'};
            case 'stop_right_earlypress' % ONLY stop-right trials when response were generated BEFORE stop signal displays
                eventtype = {'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'go_leftpress' % all left presses on go trials
                eventtype = {'Go_Left_Left';'Go_Right_Left'};
            case 'go_rightpress' % all right presses on go trials
                eventtype = {'Go_Left_Right';'Go_Right_Right'};
            case 'corr_go_left' % this flag is redundant
                eventtype = {'Go_Left_Left'};
            case 'incorr_go_left' % right response or no response on go-left trials
                eventtype = {'Go_Left_Right','Go_Left_NoResp'};
            case 'incorr_go_leftpress' % left response on go-right trials; redundant
                eventtype = {'Go_Right_Left'};
            case 'incorr_go_left_noresp' % no response on go-left trials; redundant
                eventtype = {'Go_Left_NoResp'};
            case 'corr_go_right' % redundant
                eventtype = {'Go_Right_Right'};
            case 'incorr_go_right' % left response or no response on go-right trials
                eventtype = {'Go_Right_Left','Go_Right_NoResp'};
            case 'incorr_go_rightpress' % right response on go-left trials; redundant
                eventtype = {'Go_Left_Right'};
            case 'incorr_go_right_noresp' % no response on go-right trials; redundant
                eventtype = {'Go_Right_NoResp'};
            case 'corr_stop_left' % redundant
                eventtype = {'Stop_Left_NoResp'};
            case 'incorr_stop_left' % any response on stop-left trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Left_Left','Stop_Left_Right','Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop'};
            case 'incorr_stop_left_noearlypress' % exclude stop-left trials when responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_Left','Stop_Left_Right'};
            case 'incorr_stop_left_earlypress' % ONLY stop-left trials when response were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_LeftBeforeStop','Stop_Left_RightBeforeStop'};
            case 'corr_stop_right' % redundant
                eventtype = {'Stop_Right_NoResp'};
            case 'incorr_stop_right' % any response on stop-right trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Right_Left','Stop_Right_Right','Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'incorr_stop_right_noearlypress' % exclude stop-right trials when responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Right_Left','Stop_Right_Right'};
            case 'incorr_stop_right_earlypress' % ONLY stop-right trials when response were generated BEFORE stop signal displays
                eventtype = {'Stop_Right_LeftBeforeStop','Stop_Right_RightBeforeStop'};
            case 'stop_leftpress' % any left responses on stop trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Left_Left','Stop_Right_Left','Stop_Left_LeftBeforeStop','Stop_Right_LeftBeforeStop'};
            case 'stop_leftpress_noearlypress' % exclude when left responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_Left','Stop_Right_Left'};
            case 'stop_leftpress_earlypress' % include ONLY when left responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_LeftBeforeStop','Stop_Right_LeftBeforeStop'};
            case 'stop_rightpress' % any right responses on stop trials (including when responses were generated BEFORE stop signal displays)
                eventtype = {'Stop_Left_Right','Stop_Right_Right','Stop_Left_RightBeforeStop','Stop_Right_RightBeforeStop'};
            case 'stop_rightpress_noearlypress' % exclude when right responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_Right','Stop_Right_Right'};
            case 'stop_rightpress_earlypress' % include ONLY when right responses were generated BEFORE stop signal displays
                eventtype = {'Stop_Left_RightBeforeStop','Stop_Right_RightBeforeStop'};
            case 'abort' % all early press trials
                eventtype = {'Abort_Early_Left','Abort_Early_Right'};
%           case 'noresp' % does not apply since no-response trials are
%                       explicitly correct or incorrect based on trial type
            otherwise
                eventtype = {eventtype};
        end
    case 'BDMAuctionTask'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_trialfun_hcnlab_neuroomega.'])
    case 'RatingTask'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega>reformat_eventtype.'])
    case 'RL4Task'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega>reformat_eventtype.'])
    case 'FaceWordLocalizerTask'
        ft_error(['Task ' task ' not complete. Please add event type mappings ' ...
            'to ft_read_event_hcnlab_neuroomega>reformat_eventtype.'])
    otherwise
        ft_error(['Task ' task ' not recognized. Please check or add to ' ...
            'ft_read_event_hcnlab_neuroomega.'])
end
end