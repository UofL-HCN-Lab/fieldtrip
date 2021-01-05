% eeg: plotting, visualization, averaging
%%
restoredefaultpath
addpath('C:\Users\HCN_lab\Documents\MATLAB\toolboxes\fieldtrip-20190317\')
ft_defaults
addpath(genpath('C:\Users\HCN_lab\Documents\MATLAB\FieldtripBootcamp\'))

%%
addpath('C:\Users\HCN_lab\Documents\MATLAB\FieldtripBootcamp\Fieldtrip\preprocessing_erp')
cfg              = [];
cfg.trialfun     = 'trialfun_affcog';
cfg.headerfile   = fullfile('Fieldtrip\preprocessing_erp\s04.vhdr');
cfg.datafile     = fullfile('Fieldtrip\preprocessing_erp\s04.eeg');
cfg = ft_definetrial(cfg);

%%
% Baseline-correction options
cfg.demean          = 'yes';
cfg.baselinewindow  = [-0.2 0];

% Fitering options
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 100;

% Re-referencing options - see explanation above
cfg.implicitref   = 'LM';
cfg.reref         = 'yes';
cfg.refchannel    = {'LM' 'RM'};

data = ft_preprocessing(cfg);

%%
% EOGV channel
cfg              = [];
cfg.channel      = {'53' 'LEOG'};
cfg.reref        = 'yes';
cfg.implicitref  = []; % this is the default, we mention it here to be explicit
cfg.refchannel   = {'53'};
eogv             = ft_preprocessing(cfg, data);


% only keep one channel, and rename to eogv
cfg              = [];
cfg.channel      = 'LEOG';
eogv             = ft_selectdata(cfg, eogv);
eogv.label       = {'eogv'};

% EOGH channel
cfg              = [];
cfg.channel      = {'57' '25'};
cfg.reref        = 'yes';
cfg.implicitref  = []; % this is the default, we mention it here to be explicit
cfg.refchannel   = {'57'};
eogh             = ft_preprocessing(cfg, data);

% only keep one channel, and rename to eogh
cfg              = [];
cfg.channel      = '25';
eogh             = ft_selectdata(cfg, eogh);
eogh.label       = {'eogh'};

%%
% only keep all non-EOG channels
cfg         = [];
cfg.channel = setdiff(1:60, [53, 57, 25]);        % you can use either strings or numbers as selection
data        = ft_selectdata(cfg, data);

% append the EOGH and EOGV channel to the 60 selected EEG channels
cfg  = [];
data = ft_appenddata(cfg, data, eogv, eogh);

%%
cfg        = [];
cfg.method = 'channel';
ft_rejectvisual(cfg, data)

%%
cfg          = [];
cfg.method   = 'summary';
cfg.layout   = 'mpi_customized_acticap64.mat';  % for plotting individual trials
cfg.channel  = [1:60];                          % do not show EOG channels
data_clean   = ft_rejectvisual(cfg, data);

%%
% use ft_timelockanalysis to compute the ERPs
cfg = [];
cfg.trials = find(data_clean.trialinfo==1);
task1 = ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.trials = find(data_clean.trialinfo==2);
task2 = ft_timelockanalysis(cfg, data_clean);

cfg = [];
cfg.layout = 'mpi_customized_acticap64.mat';
cfg.interactive = 'yes';
cfg.showoutline = 'yes';
ft_multiplotER(cfg, task1, task2)