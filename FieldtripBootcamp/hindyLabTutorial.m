%% set paths
restoredefaultpath
addpath('/usr/local/MATLAB/toolboxes/fieldtrip/')
ft_defaults
addpath(genpath('/usr/local/MATLAB/toolboxes/hcn_toolbox'))
addpath(genpath('/usr/local/MATLAB/toolboxes/fieldtrip/FieldtripBootcamp/'))

datapath = fullfile('/home/adamg/Documents/MATLAB/Experiments/HindyLab/');
datafile = 'Meredith~ Mari_28950aeb-9083-442a-ab35-8f2874e72111.edf';
datafilepath = fullfile(datapath,datafile);

%% load all data from TRG (or select channel) without looking at events
%%%%%%%%%% THIS WILL LIKELY TAKE FOREVER SO DO NOT TRY THIS %%%%%%%%%%%%%%%
cfg = [];
cfg.dataset = datafilepath;
cfg.continuous = 'yes';
cfg.channel = 'TRIG'; % or select channel name
cfg.trl = [46080001 46080001+5*60*2048 0];

data = ft_preprocessing(cfg);



% visually inspect the data
cfg            = [];
cfg.viewmode   = 'vertical';
ft_databrowser(cfg, data);


%%
cfg = [];
cfg.dataset = datafilepath;
cfg.continuous = 'yes';
cfg.channel = 'TRIG'; % or select channel name


t=1;
c1=1;
c2 = 1;
thresh = 1;
trigsFound = [];
trigTimeStarts = [];
while (t-1)*5*60*2048+1<= 176745472
    disp(['iteration ' num2str(t) '. time = ' num2str(((t-1)*5*60*2048+1)/2048/60)]);
    if t*5*60*2048>176745472
        stop = 176745472;
    else
        stop = t*5*60*2048;
    end
    cfg.trl = [(t-1)*5*60*2048+1 stop 0];
    data = ft_preprocessing(cfg);
    
    if any(abs(diff(data.trial{1})))
        trigsFound(c2) = 1;
        trigTimeStarts(c1) = (t-1)*5*60*2048+1;
        display(['trigger found at ' num2str(trigTimeStarts(c1))]);
        c1=c1+1;
    else
        trigsFound(t) = 0;
    end
    t =t+ 1;
    c2=c2+1;
end

cfg = [];
cfg.dataset = datafilepath;
cfg.continuous = 'yes';
cfg.channel = 'all'; % or select channel name
cfg.trl = [trigTimeStarts(1) trigTimeStarts(end)+15*60*2048 0];

data = ft_preprocessing(cfg);
dataToSave = data.trial{1};
savefilename = 'compressedFil.edf';

hdr = ft_read_header(datafilepath,'headerformat','edf');
write_edf(savefilename, hdr, dataToSave);

%% let's try reading events using default method
event = ft_read_event(datafilepath,'detect flank',[],'dataformat','edf','chanindx',145);