% functional workflow: plotting, visualization, averaging
%%
restoredefaultpath
addpath('C:\Users\HCN_lab\Documents\MATLAB\toolboxes\fieldtrip-20190317\')
ft_defaults
addpath(genpath('C:\Users\HCN_lab\Documents\MATLAB\FieldtripBootcamp\'))

%%
subjID = fullfile('Zenodo','SubjectUCI29');
load([subjID '_reref.mat'], 'reref');
load([subjID '_elec_acpc_fr.mat'], 'elec_acpc_fr');
pial_lh = ft_read_headshape('freesurfer/surf/lh.pial');
pial_lh.coordsys = 'acpc';

cfg            = [];
cfg.method     = 'mtmconvol';
cfg.toi        = -.3:0.01:.8;
cfg.foi        = 5:5:200;
cfg.t_ftimwin  = ones(length(cfg.foi),1).*0.2;
cfg.taper      = 'hanning';
cfg.output     = 'pow';
cfg.keeptrials = 'no';
freq = ft_freqanalysis(cfg, reref);

cfg            = [];
cfg.headshape  = pial_lh;
cfg.projection = 'orthographic';
cfg.channel    = {'LPG*', 'LTG*'};
cfg.viewpoint  = 'left';
cfg.mask       = 'convex';
cfg.boxchannel = {'LTG30', 'LTG31'};
lay = ft_prepare_layout(cfg, freq);

cfg              = [];
cfg.baseline     = [-.3 -.1];
cfg.baselinetype = 'relchange';
freq_blc = ft_freqbaseline(cfg, freq);

cfg             = [];
cfg.layout      = lay;
cfg.showoutline = 'yes';
ft_multiplotTFR(cfg, freq_blc);

%%
cfg             = [];
cfg.frequency   = [70 150];
cfg.avgoverfreq = 'yes';
cfg.latency     = [0 0.8];
cfg.avgovertime = 'yes';
freq_sel = ft_selectdata(cfg, freq_blc);

cfg              = [];
cfg.funparameter = 'powspctrm';
cfg.funcolorlim  = [-.5 .5];
cfg.method       = 'surface';
cfg.interpmethod = 'sphere_weighteddistance';
cfg.sphereradius = 8;
cfg.camlight     = 'no';

%%
ft_sourceplot(cfg, freq_sel, pial_lh);
view([-90 20]);
material dull;
lighting gouraud;
camlight;
ft_plot_sens(elec_acpc_fr);

