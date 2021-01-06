function [dat, label, time, cfg] = ft_preprocfun_hcnlab_mua(dat, label, time, cfg, begpadding, endpadding)

% set the defaults for the signal processing options
cfg.bpfilter =             ft_getopt(cfg, 'bpfilter', 'no');
cfg.bpfiltord =            ft_getopt(cfg, 'bpfiltord', []);
cfg.bpfilttype =           ft_getopt(cfg, 'bpfilttype', 'but');
if strcmp(cfg.bpfilttype, 'firws'), cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'onepass-zerophase'); else, cfg.bpfiltdir = ft_getopt(cfg, 'bpfiltdir', 'twopass'); end
cfg.bpinstabilityfix =     ft_getopt(cfg, 'bpinstabilityfix', 'no');
cfg.bpfiltdf =             ft_getopt(cfg, 'bpfiltdf', []);
cfg.bpfiltwintype =        ft_getopt(cfg, 'bpfiltwintype', 'hamming');
cfg.bpfiltdev =            ft_getopt(cfg, 'bpfiltdev', []);
cfg.plotfiltresp =         ft_getopt(cfg, 'plotfiltresp', 'no');
cfg.usefftfilt =           ft_getopt(cfg, 'usefftfilt', 'no');
cfg.medianfilter =         ft_getopt(cfg, 'medianfilter ', 'no');
cfg.medianfiltord =        ft_getopt(cfg, 'medianfiltord', 9);
% cfg.conv =                 ft_getopt(cfg, 'conv', 'no');
% cfg.resample =             ft_getopt(cfg, 'resample', '');

cfg.demean   =             ft_getopt(cfg, 'demean', 'no');
cfg.lpfilter =             ft_getopt(cfg, 'lpfilter', 'no');
cfg.lpfiltord =            ft_getopt(cfg, 'lpfiltord', []);
cfg.lpfilttype =           ft_getopt(cfg, 'lpfilttype', 'but');
if strcmp(cfg.lpfilttype, 'firws'), cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'onepass-zerophase'); else, cfg.lpfiltdir = ft_getopt(cfg, 'lpfiltdir', 'twopass'); end
cfg.lpinstabilityfix =     ft_getopt(cfg, 'lpinstabilityfix', 'no');
cfg.lpfiltdf =             ft_getopt(cfg, 'lpfiltdf',  []);
cfg.lpfiltwintype =        ft_getopt(cfg, 'lpfiltwintype', 'hamming');
cfg.lpfiltdev =            ft_getopt(cfg, 'lpfiltdev', []);

% compute fsample
fsample = round(1./nanmean(diff(time)));

if nargin<5 || isempty(begpadding)
  begpadding = 0;
end
if nargin<6 || isempty(endpadding)
  endpadding = 0;
end

% set the defaults for the rereferencing options
cfg.reref =                ft_getopt(cfg, 'reref', 'no');
cfg.refchannel =           ft_getopt(cfg, 'refchannel', {});
cfg.refmethod =            ft_getopt(cfg, 'refmethod', 'avg');
cfg.implicitref =          ft_getopt(cfg, 'implicitref', []);


% test whether the MATLAB signal processing toolbox is available
if strcmp(cfg.medianfilter, 'yes') && ~ft_hastoolbox('signal')
  ft_error('median filtering requires the MATLAB signal processing toolbox');
end


% do a sanity check on the rereferencing/montage
if ~strcmp(cfg.reref, 'no') && ~strcmp(cfg.montage, 'no')
  ft_error('cfg.reref and cfg.montage are mutually exclusive')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do the rereferencing in case of EEG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(cfg.implicitref) && ~any(match_str(cfg.implicitref,label))
  label = {label{:} cfg.implicitref};
  dat(end+1,:) = 0;
end

if strcmp(cfg.reref, 'yes')
  if strcmp(cfg.refmethod, 'bipolar')
    % this is implemented as a montage that the user does not get to see
    tmpcfg = keepfields(cfg, {'refmethod', 'implicitref', 'refchannel', 'channel'});
    tmpcfg.showcallinfo = 'no';
    montage = ft_prepare_montage(tmpcfg);
    % convert the data temporarily to a raw structure
    tmpdata.trial = {dat};
    tmpdata.time  = {time};
    tmpdata.label = label;
    % apply the montage to the data
    tmpdata = ft_apply_montage(tmpdata, montage, 'feedback', 'none');
    dat   = tmpdata.trial{1}; % the number of channels can have changed
    label = tmpdata.label;    % the output channels can be different than the input channels
    clear tmpdata
  else
    % mean or median based derivation of specified or all channels
    cfg.refchannel = ft_channelselection(cfg.refchannel, label);
    refindx = match_str(label, cfg.refchannel);
    if isempty(refindx)
      ft_error('reference channel was not found')
    end
    dat = ft_preproc_rereference(dat, refindx, cfg.refmethod);
  end
end

if ~strcmp(cfg.montage, 'no') && ~isempty(cfg.montage)
  % convert the data temporarily to a raw structure
  tmpdata.trial = {dat};
  tmpdata.time  = {time};
  tmpdata.label = label;
  % apply the montage to the data
  tmpdata = ft_apply_montage(tmpdata, cfg.montage, 'feedback', 'none');
  dat   = tmpdata.trial{1}; % the number of channels can have changed
  label = tmpdata.label;    % the output channels can be different than the input channels
  clear tmpdata
end


if any(any(isnan(dat)))
    % filtering is not possible for at least a selection of the data
    ft_warning('data contains NaNs, no filtering or preprocessing applied');
    
else
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do the filtering on the padded data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    transforms = cfg.transforms;
    for t=1:length(transforms)
        trans = transforms{t};
        if contains(trans,'bandpass')
            cfg.bpfilter = 'yes';
            cfg.bpfreq = [300 min([3000 fsample])];
            idx = strfind(trans,'#');
            if isempty(idx), order = 1; else, order = str2double(trans(idx+1:end)); end
            cfg.bpfiltord = order;
            dat = ft_preproc_bandpassfilter(dat, fsample, cfg.bpfreq, cfg.bpfiltord, cfg.bpfilttype, cfg.bpfiltdir, cfg.bpinstabilityfix, cfg.bpfiltdf, cfg.bpfiltwintype, cfg.bpfiltdev, cfg.plotfiltresp, cfg.usefftfilt);
        elseif contains(trans,'clip')
            idx = strfind(trans,'#');
            thresh = str2double(trans(idx+1:end));
            meanSig = mean(dat,2);
            sdSig = std(dat,0,2);
            
            pos = repmat(meanSig+thresh*sdSig,[1 size(dat,2)]);
            neg = repmat(meanSig-thresh*sdSig,[1 size(dat,2)]);
            idx1 = dat>pos;
            idx2 = dat<neg;
            dat(idx1) = pos(idx1);
            dat(idx2) = neg(idx2);
        elseif contains(trans,'volthresh')
            idx = strfind(trans,'#');
            thresh = str2double(trans(idx+1:end));
            if thresh<0
                dat = double(dat<thresh);
            else
                dat = double(dat>thresh);
            end
        elseif contains(trans,'sd') & contains(trans,'thresh')
            idx = strfind(trans,'#');
            thresh = str2double(trans(idx+1:end));
            meanSig = mean(dat,2);
            sdSig = std(dat,0,2);
            if contains(trans,'pos')
                dat = double(dat>repmat(meanSig+thresh*sdSig,[1 size(dat,2)]));
            elseif contains(trans,'neg')
                dat = double(dat<repmat(meanSig-thresh*sdSig,[1 size(dat,2)]));
            else
                idx1 = dat>repmat(meanSig+thresh*sdSig,[1 size(dat,2)]);
                idx2 = dat<repmat(meanSig-thresh*sdSig,[1 size(dat,2)]);
                dat = double(idx1 || idx2);
            end
        elseif contains(trans,'gauss')
            idx = strfind(trans,'#');
            sdms = str2double(trans(idx+1:end));
            sdsamp = fix(sdms/1000*fsample);
            x = -5*sdsamp:1:5*sdsamp;
            pd = makedist('Normal','mu',0,'sigma',sdsamp);
            gaussy = pdf(pd,x);
            gaussy = (gaussy/sum(gaussy));
            if ~rem(length(gaussy),2)
                gaussy = [gaussy 0];
            end
            dat = convn(dat, gaussy, 'same');
        elseif contains(trans,'square')
            dat = double(dat.^2);
        elseif contains(trans,'low')
            idx = strfind(trans,'#');
            idx2 = strfind(trans,'order');
            cfg.lpfreq = str2double(trans(idx(1)+1:idx2-1));
            cfg.lpfiltord = str2double(trans(idx(2)+1:end));
            dat = ft_preproc_lowpassfilter(dat, fsample, cfg.lpfreq, cfg.lpfiltord, cfg.lpfilttype, cfg.lpfiltdir, cfg.lpinstabilityfix, cfg.lpfiltdf, cfg.lpfiltwintype, cfg.lpfiltdev, cfg.plotfiltresp, cfg.usefftfilt);
        elseif contains(trans,'downsample')
            idx = strfind(trans,'#');
            resampleFs = str2double(trans(idx+1:end));
            begpadding = round(begpadding/fsample*resampleFs);
            endpadding = round(endpadding/fsample*resampleFs);
            [dat             ] = ft_preproc_resample(dat,  fsample, resampleFs, 'resample');
            [time, ~, fsample] = ft_preproc_resample(time, fsample, resampleFs, 'resample');
        else % contains(trans,'sqrt')
            dat(dat<0) = 0;
            dat = sqrt(double(dat));
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove the filter padding and do the preprocessing on the remaining trial data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if begpadding~=0 || endpadding~=0
  dat = ft_preproc_padding(dat, 'remove', begpadding, endpadding);
  if strcmp(cfg.demean, 'yes') || nargout>2
    time = ft_preproc_padding(time, 'remove', begpadding, endpadding);
  end
end
end