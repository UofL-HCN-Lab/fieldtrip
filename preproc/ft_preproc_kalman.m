function [data, cfg] = ft_preproc_kalman(data,cfg)

% reverse data so artifacts have positive peaks
data = -data;

% set the defaults
cfg.interactive         = ft_getopt(cfg, 'interactive', 'false');
cfg.arorders            = ft_getopt(cfg, 'arorders', [2,3,4,7,9]);
cfg.bestarmdlix              = ft_getopt(cfg, 'bestarmdlix', []);
cfg.bestarmodel              = ft_getopt(cfg, 'bestarmodel', []);
cfg.oeorders            = ft_getopt(cfg, 'oeorders', [2 2 1;
    3 3 1;
    3 4 1;
    5 5 1]);
cfg.arftct_zeros        = ft_getopt(cfg, 'arftct_zeros', []);
cfg.bestoemodel              = ft_getopt(cfg, 'bestoemodel', []);
cfg.bestoemdlix              = ft_getopt(cfg, 'bestoemdlix', []);
cfg.frst_artfct_endix        = ft_getopt(cfg, 'frst_artfct_endix', []);
cfg.varts               = ft_getopt(cfg, 'varts', 2.5*logspace(-1,3,5));
cfg.alphas              = ft_getopt(cfg, 'alphas', 5*logspace(-3,1,5));
cfg.bestvart               = ft_getopt(cfg, 'bestvart', []);
cfg.bestalpha              = ft_getopt(cfg, 'bestalpha',[]);
cfg.bestarmdlix = ft_getopt(cfg, 'bestarmdlix', []);
cfg.arftct_zeros = ft_getopt(cfg, 'arftct_zeros', []);
cfg.bestoemdlix = ft_getopt(cfg, 'bestoemdlix', []);
cfg.frst_artfct_endix = ft_getopt(cfg, 'frst_artfct_endix', []);

% check if 1st bad channel was already preprocessed for this trial (i.e.,
% all kalman params are non-empty)
params_found_prev = [
    ~isempty(cfg.bestarmdlix);
    ~isempty(cfg.arftct_zeros);
    ~isempty(cfg.bestoemdlix);
    ~isempty(cfg.frst_artfct_endix);
    ~isempty(cfg.bestvart);
    ~isempty(cfg.bestalpha);
    ];

if ~all(params_found_prev)
    % use previously determined parameters (if available in notes.txt)
    resp = input('Do you have previously determined parameters for this trial (yes/no)?\n','s');
    
    if strncmpi('yes',resp,length(resp))
        if isempty(cfg.bestarmdlix)
            cfg.bestarmdlix = input('Input best AR model index. Enter empty vector ([]) if you do not know.\n');
        end
        if isempty(cfg.arftct_zeros)
            cfg.arftct_zeros = input('Input zero crossing indices prior to each artifact peak. \nEnter empty vector ([]) if you do not know.\n');
        end
        if isempty(cfg.bestoemdlix)
            cfg.bestoemdlix = input('Input best OE model index. Enter empty vector ([]) if you do not know.\n');
        end
        if isempty(cfg.frst_artfct_endix)
            cfg.frst_artfct_endix = input('Input sample index of the end of the first artifact peak. \nEnter empty vector ([]) if you do not know.\n');
        end
        if isempty(cfg.bestvart)
            cfg.bestvart = input('Input best var_T. Enter empty vector ([]) if you do not know.\n');
        end
        if isempty(cfg.bestalpha)
            cfg.bestalpha = input('Input best alpha. Enter empty vector ([]) if you do not know.\n');
        end
    end
end
% do AR modeling
[cfg.bestarmodel, cfg.bestarmdlix] = best_ar_model(data,cfg);

% do PE modeling
[cfg.bestoemodel, cfg.bestoemdlix, cfg.artifact_starts] = best_oe_model(data,cfg);

% get best filter parameters & clean the data
[data, cfg.bestvart, cfg.bestalpha] = find_kalman_params_and_filter(data,cfg);

% reverse data back to original
data = -data;
close all
end

function [best_mdl, bestarmdlix] = best_ar_model(data,cfg)
Nsamp = fix(cfg.padding*cfg.Fs) - 30;
time_sample = (1:Nsamp)/cfg.Fs;
ar_orders = cfg.arorders;
data_clean_sample = data(1:Nsamp)';

if ~isempty(cfg.bestarmdlix)
    bestarmdlix = cfg.bestarmdlix;
    mdl = arima('arlags',1:ar_orders(bestarmdlix),'constant',1);
    preidx = 1:mdl.P;
    estidx = (mdl.P + 1):Nsamp;
    [best_mdl,~,~,~] = estimate(mdl,data_clean_sample(estidx),'y0',data_clean_sample(preidx));
else
    ar_mdls = cell(length(ar_orders)+1,5);
    ar_mdls(1,:) = {'Model', 'Info', 'Residual', 'Estimate', 'Fit Percent'};
    
    f1 = figure('name','AR');
    f2 = figure('name','AR Zoomed');
    f3 = figure('name','Residuals');
    for m=1:length(ar_orders)
        mdl = arima('arlags',1:ar_orders(m),'constant',1);
        preidx = 1:mdl.P;
        estidx = (mdl.P + 1):Nsamp;
        [mdl_fit,~,~,info] = estimate(mdl,data_clean_sample(estidx),'y0',data_clean_sample(preidx));
        mdl_resid = infer(mdl_fit,data_clean_sample(estidx),'Y0',data_clean_sample(preidx));
        mdl_yhat = data_clean_sample(estidx) - mdl_resid;
        mdl_prct = 1-goodnessOfFit(mdl_yhat,data_clean_sample(estidx),'nrmse');
        
        ar_mdls(m+1,:) = {mdl_fit, info, mdl_resid, mdl_yhat, mdl_prct};
        
        if cfg.interactive
            figure(f1);
            subplot(length(ar_orders),1,m);
            plot(time_sample(estidx),data_clean_sample(estidx),'b',time_sample(estidx),mdl_yhat,'r--','LineWidth',2);
            
            figure(f2);
            subplot(length(ar_orders),1,m);
            plot(time_sample(end-500:end),data_clean_sample(end-500:end),'b',time_sample(end-500:end),mdl_yhat(end-500:end),'r--','LineWidth',2);
            
            figure(f3);
            subplot(length(ar_orders),1,m);
            plot(mdl_yhat,mdl_resid,'.');
            
            if m==3
                figure(f1);ylabel('\muV');
                figure(f2);ylabel('\muV');
                figure(f3);ylabel('\muV');
            end
        end
    end
    if cfg.interactive
        figure(f1);xlabel('Time (s)');
        figure(f2);ylabel('Time (s)');
        figure(f3);ylabel('\muV');
        
        f4=figure();
        plot([ar_mdls{2:end,5}],'o-');
        ylabel('Fit');
        xlabel('AR Model');
        
        bestarmdlix = input('input index for best AR model.\n');
        
        close([f1 f2 f3 f4]);
    else
        error('Section untested');
        bestarmdlix = get_best_model_ix_auto(ar_mdls,5);
    end
    best_mdl = ar_mdls{1+bestarmdlix,1};
end
end

function [best_mdl, bestoemdlix, art_strts] = best_oe_model(data,cfg)
% a. 5 samples before & 35 samples after the artifact (1024 Hz)
oe_orders = cfg.oeorders;
oe_mdls = cell(size(oe_orders,1)+1,3);
oe_mdls(1,:) = {'Model', 'Estimate', 'Fit Percent'};

strt = 10;
nd = 35;
artBeg = fix(cfg.padding*cfg.Fs) - strt;
artEnd = length(data)-fix(cfg.padding*cfg.Fs) + nd;
artifact = data(artBeg:artEnd)';

% data_obj = iddata(artifact(2:end),artifact(1:end-1),1);
% u_t = zeros(size(artifact));  u_t(strt) = 1; %% THIS DOES NOT WORK


mn = mean(artifact);
sd = std(artifact);
if ~isempty(cfg.arftct_zeros)
    art_strts = cfg.arftct_zeros;
else
    if cfg.interactive
        figure;
        plot(artifact); hold on
        plot(1:length(artifact),(mn+sd)*ones(size(artifact)),'k--');
        art_strts = input('input indices of zero crossings prior to each artifact peak.\n');
        close(gcf);
    else
        idx_sd = artifact>=mn+sd;
        idx_0 = artifact>=0;
        idx = idx_0 & idx_sd;
        art_strts = find(diff(idx)==1)-1;
    end
end

% here, u_t must be aligned with artifact, so can use what is measured in
% lines 116-125 directly
u_t = zeros(size(artifact));  u_t(art_strts) = 1; %% THIS WORKS BEST

% this is for re-defining u_t in get_kalman_filt_params & clean_data, which
% uses the full vector data, so we need to re-align art_starts relative to
% the first sample in data
art_strts = art_strts + artBeg - 1;

% u_t(strt:(end-nd)) = 1/(artEnd-nd-(artBeg+start)+1); %% THIS DOES NOT WORK
data_obj = iddata(artifact,u_t,1);
if ~isempty(cfg.bestoemdlix)
    bestoemdlix = cfg.bestoemdlix;
    best_mdl = oe(data_obj,oe_orders(bestoemdlix,:));
else
    oefigs = zeros(size(oe_orders,1),1);
    for m=1:size(oe_orders,1)
        sys = oe(data_obj,oe_orders(m,:));
        [yhat,fit] = compare(data_obj,sys);
        oefigs(m) = figure;
        subplot(1,3,1);
        compare(data_obj,sys);
        subplot(1,3,2);
        foo = get(yhat);
        plot(1:length(artifact),artifact,'-b',1:length(artifact),artifact-foo.OutputData{1},'--r');
        subplot(1,3,3);
        pzmap(sys);
        
        oe_mdls(m+1,:) = {sys, yhat, fit/100};
    end
    
    if cfg.interactive
        figure();
        plot([oe_mdls{2:end,3}],'o-'); axis tight;
        ylims = get(gca,'ylim');
        set(gca,'ylim',ylims+[-0.02 0.02]);
        ylabel('Fit Percent');
        xlabel('OE Model');
        
        bestoemdlix = input('input index of best oe model.\n');
        close(gcf);
    else
        error('Section untested');
        bestoemdlix = get_best_model_ix_auto(oe_mdls,3);
    end
    close(oefigs);
    best_mdl = oe_mdls{1+bestoemdlix,1};
end
end

function best_ix = get_best_model_ix_auto(mdl_cell,fit_perc_ix)
fit_diffs_per_max = [nan diff([mdl_cell{2:end,fit_perc_ix}])];
fit_diffs_per_max = fit_diffs_per_max/fit_diffs_per_max(2);
fit_perc_max = [mdl_cell{2:end,fit_perc_ix}]/mdl_cell{end,fit_perc_ix};
idx1 = fit_diffs_per_max>=0.65*fit_diffs_per_max(2);
idx2 = fit_perc_max>=0.95*fit_perc_max(end);
idx = idx1 & idx2;
best_ix = find(idx,1,'first');
end

function [best_eeg_hat, best_var_T, best_alpha] = find_kalman_params_and_filter(data,cfg)

var_Ts = cfg.varts;
alphas = cfg.alphas;
tau = 30;
tau_0 = tau+1;
tau_not0 = setdiff(1:(2*tau+1),tau_0);
strt = 10;
nd = 35;
artBeg = fix(cfg.padding*cfg.Fs) - strt;
artEnd = length(data)-fix(cfg.padding*cfg.Fs) + nd;

tmpcfg = [];
t_s = fix(cfg.padding*cfg.Fs) + 1;
tmpcfg(1).t_s = t_s;

if ~isempty(cfg.frst_artfct_endix)
    pk = cfg.frst_artfct_endix;
    pk = pk + artBeg - 1;
    d = pk - t_s;
else
    if cfg.interactive
        figure();
        plot(data(artBeg:artEnd));
        pk = input('input sample index of the end of the first artifact peak.\n');
        close(gcf);
        pk = pk + artBeg - 1;
        d = pk - t_s;
    else
        d = fix(4/1024*cfg.Fs);
    end
end
tmpcfg.d = d;
d_tot = artEnd - t_s;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% state-space model for EEG
% x_E(t+1) = A_E*x_E(t) + G_E*e_E(t), where e_E(t)->WN(0,var_E)
% eeg(t) = C_E*x_E(t)
ar_mdl = cfg.bestarmodel;
ord = length(ar_mdl.AR);
tmpcfg.ord = ord;

A_E = [zeros(1,ord);[eye(ord-1) zeros(ord-1,1)]];
A_E(1,:) = -cell2mat(ar_mdl.AR);

G_E = zeros(1,ord); G_E(1) = 1; % need to transpose at end
C_E = zeros(1,ord); C_E(end) = 1;
tmpcfg.C_E = C_E;

var_E = ar_mdl.Variance;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% state-space model for TMS
% x_T(t+1) = A_T*x_T(t) + B_T*u(t), where u(t)-> impulse
% tms(t) = C_T*x_T(t) + nu(t), nu(t)-> WN(0,var_nu)
oe_mdl = cfg.bestoemodel;
deg = max(length(oe_mdl.B)-1,length(oe_mdl.F)-1);
tmpcfg.deg = deg;

st_sp = ss(oe_mdl);
A_T = st_sp.A;

B_T = st_sp.B'; % will be transposed back at end
C_T = st_sp.C;
tmpcfg.C_T = C_T;

var_nu = oe_mdl.NoiseVariance;
tmpcfg.var_nu = var_nu;

u_t = zeros(size(data));  u_t(cfg.artifact_starts) = 1;
tmpcfg.u_t = u_t;
%%%%%%%%
%%%%%%%% NEXT: try implementing state-space model with MATLAB's built-in
%%%%%%%% kalman function
%%%%%%%%

% uy = impulse(st_sp);
% u_t = zeros(size(data));
% u_t(t_s:(t_s+length(uy)-1)) = uy;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% setting and collecting state-space systems
RE = size(A_E,1); RT = size(A_T,1);
RT0 = RT; RE0 = RE;

A = [A_E zeros(RE0,RT0);
    zeros(RT0,RE0) A_T];
tmpcfg.A = A;

B = [zeros(1,RE0) B_T]'; % needs to be column vector
tmpcfg.B = B;
G = [G_E zeros(1,RT0)]'; % needs to be column vector

C = [C_E C_T]; % should be row vector
tmpcfg.C = C;

G_M = [G [zeros(RE0,RT0); eye(RT0)]];
tmpcfg.G_M = G_M;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% test filter over var_Ts and alphas
Nsamp = fix(cfg.padding*cfg.Fs) - 30;
eeg_hat = repmat(mean(data),[1 length(data)]);

prev_white = -inf;
prev_ind = -inf;
prev_eeg_rmse = sqrt(mean((eeg_hat(1:Nsamp) - data(1:Nsamp)).^2));
prev_pow_rmse = inf;
prev = [prev_white prev_ind 1/prev_eeg_rmse 1/prev_pow_rmse];

Q_time1 = t_s:(t_s+d);
tmpcfg.Q_time1 = Q_time1;
% % make Q_time1s go to ds_tot
% Q_time1 = ts_s:(ts_s+d_tot);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% ------------------------------ START HERE -------------------------------
% add option of no eta_time2 (no decay of artifact), fix eta_t, R_t and Rix
% in kalman_local when tmpcfg.include_nu_decay = false
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
tmpcfg.include_nu_decay = cfg.include_nu_decay;
if tmpcfg.include_nu_decay
    eta_time1 = t_s:(t_s+d);
    eta_time2 = (t_s+d+1):(t_s+d_tot);
else
    eta_time1 = [];
    eta_time2 = t_s:(t_s+d_tot) + d + 1;
end
tmpcfg.eta_time1 = eta_time1;
tmpcfg.eta_time2 = eta_time2;

best_var_T = nan;
best_alpha = nan;

% for interactive, save some data by alpha level to plot summary of all
% alpha values for each var_T
eeg_hat_byalpha = zeros(length(alphas),size(eeg_hat,2));
R_t2_byalpha = zeros(length(alphas),length(eta_time2));
artfct_byalpha = zeros(length(alphas),length(data));
Rhat_eps_CI_byalpha = zeros(length(alphas),2*tau+1);
Rhat_eps_byalpha = zeros(length(alphas),2*tau+1);
Rhat_eps_eeg_CI_byalpha = zeros(length(alphas),2*tau+1);
Rhat_eps_eeg_byalpha = zeros(length(alphas),2*tau+1);
power_hat_y_byalpha = zeros(length(alphas),floor(Nsamp/2));

eeg_hat_prevbest = nan(size(eeg_hat));
R_T2_prevbest = nan(size(eta_time2));
artfct_prevbest = nan(size(data));
Rhat_eps_CI_prevbest = nan(1,2*tau+1);
Rhat_eps_prevbest = nan(1,2*tau+1);
Rhat_eps_eeg_CI_prevbest = nan(1,2*tau+1);
Rhat_eps_eeg_prevbest = nan(1,2*tau+1);
power_hat_y_prevbest = nan(1,floor(Nsamp/2));
f_prev = figure();
f_cur1 = figure();
f_cur2 = figure();
for v=1:length(var_Ts)
    if ~isempty(cfg.bestvart) && v==1
        var_T = cfg.bestvart;
    elseif isempty(cfg.bestvart)
        var_T = var_Ts(v);
    else
        continue
    end
    
    % complete covariance matrix as a function of t relative to t_s and d
    % t_s is instant in which arifact is generated
    Q_t{1} = diag([var_E var_T*ones(1,deg)]);  % for t_s <= t <= t_s+d
    Q_t{2} = diag([var_E zeros(1,deg)]);       % otherwise
    tmpcfg.Q_t = Q_t;
    for a=1:length(alphas)
        if ~isempty(cfg.bestalpha) && a==1
            tmpcfg.alpha = cfg.bestalpha;
        elseif isempty(cfg.bestalpha)
            tmpcfg.alpha = alphas(a);
        else
            continue
        end
        [results] = kalman_local(data,tmpcfg);
        
        eeg_hat = results.eeg_hat;
        eeg_hat_byalpha(a,:) = eeg_hat;
        artfct_byalpha(a,:) = results.artfct;
        R_t = results.R_t;
        R_t2_byalpha(a,:) = R_t{2};
        
        % whiteness & independence test calculations
        Ns = nan(size(tau));
        eeg_hat_sample = eeg_hat(1,1:Nsamp);
        % whiteness data setup
        % Rhat_eps = mean(eps(t)*eps(t-tau)), where eps = eeg - eeg_hat,
        % where eeg = output of AR model for clean eeg data
        data_clean_sample = data(1,1:Nsamp);
        preidx = 1:ar_mdl.P;
        estidx = (ar_mdl.P + 1):Nsamp;
        eeg_resid = infer(ar_mdl,data_clean_sample(estidx)','Y0',data_clean_sample(preidx)')';
        eeg_sample = data_clean_sample(estidx) - eeg_resid;
        eps_sample = eeg_sample - eeg_hat_sample(estidx);
        N = length(eps_sample);
        
        
        % independence & whiteness calculations
        ix = 1;
        for t=-tau:tau
            Ns(ix) = N - abs(t);
            ix = ix+1;
        end
        Rhat_eps = xcorr(eps_sample,tau,'biased');
        Rhat_eps0 = Rhat_eps(tau_0);
        Rhat_eps_CI = 10*2.58./sqrt(N);
%         Rhat_eps_CI = 10*2.58./sqrt(Ns);% not sure Morbidi et al. (2007) came up with their threshold
                                        % in Figure 8a, seems to be off by a factor of 10 based on the 
                                        % numbers they provide, I did 10x here to be consistent with 
                                        % what it looks like in their fig
        Rhat_eps_CI_byalpha(a,:) = Rhat_eps_CI;
        Rhat_eps_byalpha(a,:) = Rhat_eps/Rhat_eps0;
        
        Rhat_eeg = xcorr(eeg_hat_sample(estidx),tau,'biased');
        Rhat_eps_eeg = xcorr(eps_sample,eeg_hat_sample(estidx),tau,'biased');
        S = sum(Rhat_eps.*Rhat_eeg);
        Rhat_eps_eeg_CI = 10*2.58*sqrt(S/N);
%         Rhat_eps_eeg_CI = 10*2.58*sqrt(S./Ns); % not sure Morbidi et al. (2007) came up with their threshold
                                        % in Figure 8b (depends on their
                                        % values of S). for consistency
                                        % with what I did for Rhat_eps_CI,
                                        % I also did 10x here...
        Rhat_eps_eeg_CI_byalpha(a,:) = Rhat_eps_eeg_CI;
        Rhat_eps_eeg_byalpha(a,:) = Rhat_eps_eeg;
        
        cur_white = sum(Rhat_eps(tau_not0)/Rhat_eps0>=-Rhat_eps_CI & Rhat_eps(tau_not0)/Rhat_eps0<= Rhat_eps_CI);
%         cur_white = sum(Rhat_eps(tau_not0)/Rhat_eps0>=-Rhat_eps_CI(tau_not0) & Rhat_eps(tau_not0)/Rhat_eps0<= Rhat_eps_CI(tau_not0));
        cur_ind = sum(Rhat_eps_eeg>=-Rhat_eps_eeg_CI & Rhat_eps_eeg<= Rhat_eps_eeg_CI);
        
        % spectal comparison of AR model of clean LFP and cur_clean
        cur_eeg_rmse = sqrt(mean((eeg_hat_sample - data_clean_sample).^2));
        
        eeg_y = fft(hann(length(data_clean_sample))'.*data_clean_sample);
        eeg_hat_y = fft(hann(length(eeg_hat_sample))'.*eeg_hat_sample);
        
        f = (0:(length(eeg_hat_sample)-1))*cfg.Fs/length(eeg_hat_sample);
        power_y = abs(eeg_y).^2/length(data_clean_sample);
        power_hat_y = abs(eeg_hat_y).^2/length(eeg_hat_sample);
        
        f = f(1:floor(length(eeg_hat_sample)/2));
        power_y = power_y(1:floor(length(eeg_hat_sample)/2));
        power_hat_y = power_hat_y(1:floor(length(eeg_hat_sample)/2));
        
        power_hat_y_byalpha(a,:) = power_hat_y;
        
        cur_pow_rmse = sqrt(mean((power_hat_y - power_y).^2));
        
        cur = [cur_white cur_ind 1/cur_eeg_rmse 1/cur_pow_rmse];
        
        % if current results better than previous results, replace
        % previous results with current
        if isempty(cfg.bestvart) || isempty(cfg.bestalpha)
            if ~cfg.interactive
                chk = prev<cur;
                if sum(chk)>=3 % at least 3 tests using current params are better than previous
                    best_eeg_hat = eeg_hat;
                    prev = cur;
                    best_var_T = var_T;
                    best_alpha = tmpcfg.alpha;
                end
            end
        end
    end
    if isempty(cfg.bestvart) || isempty(cfg.bestalpha)
        if cfg.interactive
            ix = (artBeg-40):(artEnd+20);
            % plot results from previously best found alpha & var_T
            figure(f_prev);
            subplot(5,1,1);
            plot(ix,data(ix),'b-',ix,eeg_hat_prevbest(ix),'r--');
            axis tight;
            title('Data comparison');
            subplot(5,1,2);
            plot(ix,artfct_prevbest(ix));
            axis tight;
            xlabel('Time (s)');
            ylabel('\muV');
            title('Estimated artifact');
            subplot(5,1,3);
            plot(-tau:tau,-Rhat_eps_CI_prevbest,'k--',-tau:tau,Rhat_eps_prevbest,'b',-tau:tau,Rhat_eps_CI_prevbest,'k--');
            axis tight;
            xlabel('\tau');
            ylabel('$$\hat{R}_{\epsilon}$$','interpreter','latex');
            title('Whiteness Test');
            subplot(5,1,4);
            plot(-tau:tau,-Rhat_eps_eeg_CI_prevbest,'k--',-tau:tau,Rhat_eps_eeg_prevbest,'b',-tau:tau,Rhat_eps_eeg_CI_prevbest,'k--');
            axis tight;
            xlabel('\tau');
            ylabel('$$\hat{R}_{\xi,\hat{eeg}}$$','interpreter','latex');
            title('Independence Test');
            subplot(5,1,5);
            loglog(f,power_y,'b',f,power_hat_y_prevbest,'r--');
            axis tight;
            xlabel('Frequency (Hz)');
            ylabel('Power');
            title('Power comparison');
            
            % plot comparison of cur_clean and trial_data
            figure(f_cur1);
            for a=1:length(alphas)
                subplot(length(alphas),1,a)
                plot(ix,eeg_hat_prevbest(ix),'k:','linewidth',2);
                hold on;
                plot(ix,data(ix),'b-',ix,eeg_hat_byalpha(a,ix),'r--');%,...
                %                 [eta_time1 eta_time2],[R_t{1}*ones(size(eta_time1)) R_t2_byalpha(a,:)],'k--');
                axis tight;
            end
            xlabel('Time (s)');
            ylabel('\muV');
            
            figure(f_cur2);
            plot(ix,artfct_prevbest(:,ix),'k:','linewidth',2);
            hold on;
            plot(ix,artfct_byalpha(:,ix));
            axis tight;
            legend();
            xlabel('Time (s)');
            ylabel('\muV');
            
            rm1 = figure('name','Whiteness & Independence Tests');
            for a=1:length(alphas)
                subplot(length(alphas),2,2*(a-1)+1);
                plot(-tau:tau,-Rhat_eps_CI_byalpha(a,:),'k--',-tau:tau,Rhat_eps_byalpha(a,:),'b',-tau:tau,Rhat_eps_CI_byalpha(a,:),'k--');
                axis tight;
                if a==1
                    title('Whiteness');
                end
                if a==length(alphas)
                    xlabel('\tau');
                    ylabel('$$\hat{R}_{\epsilon}$$','interpreter','latex');
                end
                subplot(length(alphas),2,2*(a-1)+2);
                plot(-tau:tau,-Rhat_eps_eeg_CI_byalpha(a,:),'k--',-tau:tau,Rhat_eps_eeg_byalpha(a,:),'b',-tau:tau,Rhat_eps_eeg_CI_byalpha(a,:),'k--');
                axis tight;
                if a==1
                    title('Independence');
                end
                if a==length(alphas)
                    xlabel('\tau');
                    ylabel('$$\hat{R}_{\xi,\hat{eeg}}$$','interpreter','latex');
                end
            end
            
            rm2 = figure('name','Spectrum Comparisons');
            for a=1:length(alphas)
                subplot(length(alphas),1,a);
                loglog(f,power_y,'b',f,power_hat_y_byalpha(a,:),'r--');
                axis tight;
            end
            xlabel('Frequency (Hz)');
            ylabel('Power');
            
            best_ix = input('input row number of best filter. if none of these options are better than previous best, return empty vector ([]).\n');
            if ~isempty(best_ix)
                eeg_hat_prevbest = eeg_hat_byalpha(best_ix,:);
                R_T2_prevbest = R_t2_byalpha(best_ix,:);
                artfct_prevbest = artfct_byalpha(best_ix,:);
                Rhat_eps_CI_prevbest = Rhat_eps_CI_byalpha(best_ix,:);
                Rhat_eps_prevbest = Rhat_eps_byalpha(best_ix,:);
                Rhat_eps_eeg_CI_prevbest = Rhat_eps_eeg_CI_byalpha(best_ix,:);
                Rhat_eps_eeg_prevbest = Rhat_eps_eeg_byalpha(best_ix,:);
                power_hat_y_prevbest = power_hat_y_byalpha(best_ix,:);
                
                best_var_T = var_T;
                best_alpha = alphas(best_ix);
                best_eeg_hat = eeg_hat_prevbest;
            end
            close([rm1 rm2]);
            clf(f_cur1);
            clf(f_cur2);
        end
    else
        best_var_T = cfg.bestvart;
        best_alpha = cfg.bestalpha;
        best_eeg_hat = eeg_hat;
    end
end
end

function [results] = kalman_local(data,cfg)
ord = cfg.ord;
deg = cfg.deg;
C = cfg.C;
alpha = cfg.alpha;
eta_time1 = cfg.eta_time1;
eta_time2 = cfg.eta_time2;
t_s = cfg.t_s;
d = cfg.d;
var_nu = cfg.var_nu;
C_E = cfg.C_E;
C_T = cfg.C_T;
Q_t = cfg.Q_t;
Q_time1 = cfg.Q_time1;
A = cfg.A;
B = cfg.B;
G_M = cfg.G_M;
u_t = cfg.u_t;

xi_tplus1 = zeros(size(data));
eeg_hat = nan(size(data));
artfct = nan(size(data));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% itialize filter (first prediction step)
xhat_tplus1_t = zeros(ord+deg,1);
P_tplus1_t = [eye(ord) zeros(ord,deg);
    zeros(deg,ord) 1e-6*eye(deg)];

% first correction step
xi_tplus1(1) = data(1) - C*xhat_tplus1_t;
R_tplus1 = 0; % first time point, R=0
K_tplus1 = P_tplus1_t*C'*((C*P_tplus1_t*C'+R_tplus1)^(-1));
xhat_tplus1_tplus1 = xhat_tplus1_t + K_tplus1*xi_tplus1(1);
P_tplus1_tplus1 = P_tplus1_t - K_tplus1*C*P_tplus1_t;

eta_t = -alpha*(eta_time2 - t_s - d);

% replace nu(t) with eta(t), model variance of eta(t) as:
R_t{1} = var_nu;
R_t{2} = var_nu*exp(eta_t);
R_t{3} = 0;

% - - - - - - - - - - - - - - - - - - - - - - - - - - -- - - - - - - - - - 
for tplus1=2:length(data)
    xhat_t_t = xhat_tplus1_tplus1;
    P_t_t = P_tplus1_tplus1;
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % current clean LFP estimate
    eeg_hat(tplus1-1) = C_E*xhat_t_t(1:ord);
    artfct(tplus1-1) = C_T*xhat_t_t(ord+1:end);
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % determine sample-dependent values
    %             u_t = (1/(eta_time2(end)-eta_time1(1))) * ismember(tplus1-1,[eta_time1 eta_time2]); % for t during artifact
    %             u_t = tplus1-1==t_s; % for t == artifact start
    %             u_t = cfg.Fs * tplus1-1==t_s; % for t == artifact start
    
    Q_tplus1 = Q_t{double(~ismember(tplus1,Q_time1))+1};
    Rix = [ismember(tplus1,eta_time1), ismember(tplus1,eta_time2), ~ismember(tplus1,[eta_time1 eta_time2])];
    Rix_ix = [1 find(eta_time2==tplus1)];
    R_tplus1 = R_t{Rix}(Rix_ix(end));
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % prediction step for estimate of next sample
    %             xhat_tplus1_t = A*xhat_t_t + B*u_t;
    xhat_tplus1_t = A*xhat_t_t + B*u_t(tplus1-1);
    P_tplus1_t = A*P_t_t*A' + G_M*Q_tplus1*G_M';
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % correction step for current estimate
    xi_tplus1(tplus1) = data(tplus1) - C*xhat_tplus1_t;
    K_tplus1 = P_tplus1_t*C'*((C*P_tplus1_t*C'+R_tplus1)^(-1));
    %                     K_t = P_t_tminus1*C'\(C*P_t_tminus1*C'+R);
    
    xhat_tplus1_tplus1 = xhat_tplus1_t + K_tplus1*xi_tplus1(tplus1);
    P_tplus1_tplus1 = P_tplus1_t - K_tplus1*C*P_tplus1_t;
end
% last estimate
xhat_t_t = xhat_tplus1_tplus1;
P_t_t = P_tplus1_tplus1;
eeg_hat(tplus1) = C_E*xhat_t_t(1:ord);
artfct(tplus1) = C_T*xhat_t_t(ord+1:end);

results.xhat_t_t = xhat_t_t;
results.P_t_t = P_t_t;
results.eeg_hat = eeg_hat;
results.artfct = artfct;
results.R_t = R_t;
end
