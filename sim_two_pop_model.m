%%-----------------------------------------------------------------------%%
% This script runs stochastic simulations of the two population model.
% This script produces plots comparing the empirical cumulative distribution
% of clone sizes with the one obtained numerically.
% For details on the model and numerical implementations please refer to
% the METHODS S1 document.
%%-----------------------------------------------------------------------%%
close all;
clear all; 
clc;

% provide plots directory and RUN
plot_dir = 'plots';

datas = {'Confetti','Red2Kras'};
for sets = 1:2
    data = datas{sets};
    if sets == 2
        chans = 2; chs = {'YFP','RFP'};
    else
        chans = 1; chs = {''};
    end
    for chan = 1:chans
        switch data
            case 'Confetti'
                model = 1;
                fss_in = 0.16;
                sigmas_s_in = 2*[0.124]/7;
                sigmas_p_in = 2*[0.03]/7;
                r_in = 0.5;%0.5
                q_in = 0.5;
                times = 7*[12 24 36 52 60 72];
            case 'Red2Kras'
                if chan == 1
                    model = 2;
                    %             chan = 1;
                    %         fss_in = 0.3; % 4week
                    fss_in = 0.12; % 4week
                    %         sigmas_s_in = 2*3.55/7;
                    %         sigmas_p_in = 2*0.6/7;
                    sigmas_s_in = 2*[13.62]/7;%2*2.55/7;
                    sigmas_p_in = 2*([1.5])/7;%2*0.6/7;
                    r_in = 0.5;%0.5
                    q_in = 0.5;%0.5*(1-params.sigma_s/(18*params.sigma_p));
                    times = 7*[1,2,4];

                    tau = 7*1;
                    taur = r_in;
                    tauq = q_in;
                    tausigma_s = 2*1.221/7;
                    tausigma_p = 2*(0.2425)/7;
                elseif chan ==2
                    model = 2;
                    %         chan = 2;
                    fss_in = 0.08;
                    r_in = 0.7;
                    q_in = r_in;
                    sigmas_s_in = 3.1/7/(2*r_in-1);
                    sigmas_p_in = 0.9/7/(2*q_in-1);
                    tau = 7*2;
                    taur = r_in;
                    tauq = q_in;
                    tausigma_s = 0.5/7/(2*taur-1);
                    tausigma_p = 0.01/7/(2*tauq-1);
                    times = 7*[1,2,4];
                end
        end

        nfs = 0;
        parameters = [];
        for sigma_s_in = sigmas_s_in
            for sigma_p_in = sigmas_p_in
                for fs = fss_in
                    nfs = nfs +1;
                    parameters(nfs,:) = [sigma_s_in,sigma_p_in,fs];
                    params.sigma_s = sigma_s_in;
                    params.sigma_p = sigma_p_in;
                    params.r = r_in;
                    params.q = q_in;
                    params.fs = fs;
                    load_figs = 1;
                    nclones = 1000;

                    ignore_clones_of_size = 1;

                    N_max = Inf;
                    n_realis = 1000;
                    fs_vs_t = [];
                    sp_clone_sizes_vs_t = {};

                    nt = 0;
                    sim_clones_size_all = [];
                    for tmax = times
                        params.tmax = tmax;
                        if load_figs == 1
                            switch data
                                case 'Confetti'
                                    fig = openfig([plot_dir,'/errorbar_csdist_wo_cs1_conf',num2str(tmax/7),'w.fig']);
                                    dataObjs = findobj(fig,'-property','YData');
                                    y1 = dataObjs(1).YData;
                                    x1 = dataObjs(1).XData;
                                case 'Red2Kras'
                                    if chan == 1
                                        fig = openfig([plot_dir,'/errorbar_csdist_wo_cs1_kras',num2str(tmax/7),'w_ch_YFP.fig']);

                                        dataObjs = findobj(fig,'-property','YData');
                                        y1 = dataObjs(1).YData;
                                        x1 = dataObjs(1).XData;
                                    else
                                        fig = openfig([plot_dir,'/errorbar_csdist_wo_cs1_kras',num2str(tmax/7),'w_ch_RFP.fig']);
                                        dataObjs = findobj(fig,'-property','YData');
                                        y1 = dataObjs(1).YData;
                                        x1 = dataObjs(1).XData;
                                    end
                            end
                        end
                        fig.Position = [000 100 250 250];
                        pause(1)
                        lims = axis();
                        counter = 0;
                        all_sim_clone_sizes = []; s_clone_sizes = []; p_clone_sizes = [];
                        sp_clone_sizes = [];
                        counter = counter + 1;
                        r = params.r;

                        r_orgininal = 1;
                        factor = (2*r_orgininal-1)/r_orgininal;
                        sigma_s = params.sigma_s;
                        q = params.q;
                        sigma_p = params.sigma_p;
                        fs =  params.fs;

                        delta_s = (2*r-1)*sigma_s;
                        delta_p =  (2*q-1)*sigma_p;
                        binedges = ignore_clones_of_size+0.5:1:100000.5; % define bins
                        cumNsims = zeros(n_realis,length(binedges)-1); % preallocate array
                        % run n_realis realizations
                        for n = 1:n_realis
                            % each realisation simulates nclones clones
                            if model == 1
                                clone_info = function_2_population_model(sigma_s,r,sigma_p,q,nclones,fs,tmax,N_max);
                            elseif model == 2
                                clone_info = function_2_population_model_ageing(sigma_s,r,sigma_p,q,tau,tausigma_s,taur,tausigma_p,tauq,nclones,fs,tmax,N_max);
                            end
                            % clone size = number of S cells + number of P cells in clone
                            sim_clones_size = clone_info.numF + clone_info.numS;

                            if ignore_clones_of_size >= 0
                                sim_clones_size(sim_clones_size <= ignore_clones_of_size) = [];
                            end
                            % construct cdf for every realisation
                            [~,cumNsim,bincents] = cumulativeProb(sim_clones_size,binedges);
                            cumNsims(n,:) = cumNsim;
                            % store the size of all clones (for KS-test later)
                            sp_clone_sizes = [sp_clone_sizes;clone_info.numF,clone_info.numS];

                            ind = find(clone_info.numF_init);
                            s_clone_sizes = [s_clone_sizes;clone_info.numF(ind)];
                            ind = find(clone_info.numS_init);
                            p_clone_sizes = [p_clone_sizes;clone_info.numS(ind)];

                            all_sim_clone_sizes = [all_sim_clone_sizes;sim_clones_size];
                        end
                        % plot numerical cdf (average +/- SD)
                        cdfMean = mean(cumNsims,1); cdfSD = std(cumNsims,[],1);

                        hold on
                        plot(bincents,cdfMean,'-','Color','k','LineWidth',1.5)
                        plot(bincents,cdfMean+cdfSD,'--','Color','k','LineWidth',1,'HandleVisibility','off')
                        plot(bincents,cdfMean-cdfSD,'--','Color','k','LineWidth',1)

                        set(gca,'YScale','log')
                        xlabel('Clone size (cells)'); ylabel('Cumulative probability')
                        legend({[num2str(tmax/7),' week ',data,' ',chs{chan}],'Biexp. fit (mean)','Biexp. fit (SD)'},'Location','northeast'); legend boxoff
                        set(gcf,'color','w'); box on
                        axis(lims)
                    end
                end
            end
        end
    end
end
%% Functions
function clone_info = function_2_population_model(sigma_f,r,sigma_s,q,nclones,prop_s,t_max,N_max)
% Simulations of the two population model,  METHODS S1 Eq. (1)
%
% F : fast cycling cells
% S : slow cycling cells
% processes:
% sigma_f   : F -> F+F Prob. r (renewal through symmetric division)
% sigma_f   : F -> 0   Prob 1-r (loss)
% sigma_s : S -> S+S Prob. q (renewal through symmetric division)
% sigma_s : S -> 0 Prob 1-q (loss)
% ----------------------------------------------------------------------- %
% Inputs:
% sigma_f : F cycling rate
% r       : F duplication probability
% sigma_s : S cycling rate
% q       : S duplication probability
% nclones : Number of realizations (each realization simulates a single clone)
% prop_f  : Proportion of Nreps that are initialised with a single F cell
%           i.e., S(t=0)=1, P(t=0)=0, the remaining 1-prop_f realisations are 
%           initialized with a single S cell i.e., S(t=0)=0, P(t=0)=1.
% t_max   : Maximum simulation time (equates to the chase time)
% N_max   : Upper bound for the total number of particles F+S in the system,
%           can be set to Inf.
%
% A single realisation terminates if t_max or N_max are reached, or if the
% total particle count S+P is zero.
% If S+P is zero, then the realization in re-run.
%
% Outputs:
% clone_info : table where each row shows the results for a single 
%               realisation. Columns show total simulation time, initial 
%               number of F (numF_init) and S (numS_init) cells, and final 
%               number of F (numF) and S (numS) cells. 
% -------------------------------------------------------------------------
clone_sizes = zeros(nclones,5); % initialize array to store results
rep = 0;

while rep < nclones % run until nclones clones are recorded
    % initial condition:
    if rep <= prop_s*nclones
        numF_init = 1; numS_init = 0;
    else
        numF_init = 0; numS_init = 1;
    end
    % initialize variables
    numF = numF_init;  % number of stem cells
    numS = numS_init;  % number of progenitors
    t = 0;             % initial time
    % run realisation until one of the following conditions is met: 
    % (i) reaching maximum time, (ii) reaching max number of particles, or
    % (iii) the clone dies (numS + numP == 0).
    while t < t_max && (numF + numS) < N_max && (numF + numS) > 0
        % calculate propensity functions:
        wf1 = sigma_f*r*numF;
        wf2 = sigma_f*(1-r)*numF;
        ws1 = sigma_s*q*numS;
        ws2 = sigma_s*(1-q)*numS;
        w = wf1 + wf2 + ws1 + ws2;
        % update time until the next event
        dt =  - log(1-rand()) / w;
        t = t + dt;
        % choose an action
        ran = rand();
        if ran <= wf1/w              % F -> F + F
            numF = numF + 1;
        elseif ran <= (wf1+wf2)/w    % F -> 0
            numF = numF - 1;
        elseif ran <= (wf1+wf2+ws1)/w % S -> S + S
            numS = numS + 1;
        elseif ran <= (wf1+wf2+ws2+ws2)/w % S -> 0
            numS = numS - 1;
        end  
    end
   
    if (numF + numS) >= 0 % if the clone is "proliferative"
        % store clone size, and initial condition information
        rep = rep + 1;
        clone_sizes(rep,:) = [t, numF_init, numS_init, numF, numS];
    end
end
% return table:
clone_info = array2table(clone_sizes,...
    'VariableNames',{'t','numF_init','numS_init','numF','numS'});
end % function_2_population_model

%% Functions
function clone_info = function_2_population_model_ageing(sigma_f,r,sigma_s,q,tau,tausigma_f,taur,tausigma_s,tauq,nclones,prop_f,t_max,N_max)
% Simulations of the two population model,  METHODS S1 Eq. (1), considering
% a change in the proliferation dynamic after time tau.
%
% F : fast cycling cells
% S : slow cycling cells
% processes for time t in [0,tau]:
% sigma_f   : F -> F+F Prob. r (renewal through symmetric division)
% sigma_f   : F -> 0   Prob 1-r (loss)
% sigma_s : S -> S+S Prob. q (renewal through symmetric division)
% sigma_s : S -> 0 Prob 1-q (loss)
% processes for time t in (tau,inf]:
% tausigma_f   : F -> F+F Prob. taur (renewal through symmetric division)
% tausigma_f   : F -> 0   Prob 1-taur (loss)
% tausigma_s : S -> S+S Prob. tauq (renewal through symmetric division)
% tausigma_s : S -> 0 Prob 1-tauq (loss)
% ----------------------------------------------------------------------- %
% Inputs:
% sigma_f    : F cycling rate
% r          : F duplication probability
% sigma_s    : S cycling rate
% q          : S duplication probability
% tau        : time at which rates change
% tausigma_f : new F cycling rate (after t=tau)
% taur       : new F duplication probability
% tausigma_s : new S cycling rate
% tauq       : new S duplication probability
% nclones : Number of realizations (each realization simulates a single clone)
% prop_f  : Proportion of Nreps that are initialised with a single F cell
%           i.e., S(t=0)=1, P(t=0)=0, the remaining 1-prop_f realisations are 
%           initialized with a single S cell i.e., S(t=0)=0, P(t=0)=1.
% t_max   : Maximum simulation time (equates to the chase time)
% N_max   : Upper bound for the total number of particles F+S in the system,
%           can be set to Inf.
%
% A single realisation terminates if t_max or N_max are reached, or if the
% total particle count S+P is zero.
% If S+P is zero, then the realization in re-run.
%
% Outputs:
% clone_info : table where each row shows the results for a single 
%               realisation. Columns show total simulation time, initial 
%               number of F (numF_init) and S (numS_init) cells, and final 
%               number of F (numF) and S (numS) cells. 
% -------------------------------------------------------------------------
clone_sizes = zeros(nclones,5); % initialize array to store results
rep = 0;
while rep < nclones % run until Nreps clones are recorded
    % initial condition:
    if rep <= prop_f*nclones
        numF_init = 1; numS_init = 0;
    else
        numF_init = 0; numS_init = 1;
    end
    % initialize variable
    numF = numF_init;  % number of stem cells
    numS = numS_init;  % number of progenitors
    t = 0;             % initial time
    % run realisation until reaching maximum time or number of particles.
    % or the clone dies (numS + numP == 0).
    while t < t_max && (numF + numS) < N_max && (numF + numS) > 0
        % calculate propensity functions:
        if t < tau
            ws1 = sigma_f*r*numF;
            ws2 = sigma_f*(1-r)*numF;
            wp1 = sigma_s*q*numS;
            wp2 = sigma_s*(1-q)*numS;
        else
            ws1 = tausigma_f*taur*numF;
            ws2 = tausigma_f*(1-taur)*numF;
            wp1 = tausigma_s*tauq*numS;
            wp2 = tausigma_s*(1-tauq)*numS;
        end
        w = ws1 + ws2 + wp1 + wp2;
        % update time until the next event
        if w>0
        dt =  - log(1-rand()) / w;
        t = t + dt;
        % choose an action
        ran = rand();
        if ran <= ws1/w              % F -> F + F
            numF = numF + 1;
        elseif ran <= (ws1+ws2)/w    % F -> 0
            numF = numF - 1;
        elseif ran <= (ws1+ws2+wp1)/w % S -> S + S
            numS = numS + 1;
        elseif ran <= (ws1+ws2+wp2+wp2)/w % S -> 0
            numS = numS - 1;
        end  
        end
    end
   
    if (numF + numS) >= 0 % if the clone is "proliferative"
        % store clone size, and initial condition information
        rep = rep + 1;
        clone_sizes(rep,:) = [t, numF_init, numS_init, numF, numS];
    end
end
% return table:
clone_info = array2table(clone_sizes,...
    'VariableNames',{'t','numF_init','numS_init','numF','numS'});
end % function_2_population_model_ageing

function [N,cumN,bincents] = cumulativeProb(data,binedges)
% Returns the probabilities N and cumulative probabilities cumN, and
% corresponding bincentres for the Nx1 data array. 
if nargin == 1
binedges = 0.5:1:10^5+0.5; 
end
bincents = (binedges(2:end) + binedges(1:end - 1))/2;
N = histcounts(data,binedges,'Normalization','pdf');
cumN = 1 - cumsum(N)./sum(N);
end
