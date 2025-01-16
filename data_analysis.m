close all;
clear all;
clc;
if ispc; slsh = '\'; else; slsh = '/'; end
root_path = '';
%
% This script load the data in datadir and plots the panels indicated at
% the begining of each section.
% The scrpt is meant to be run top to bottom, as some panels depend on some
% preprocesing of the data that is not in the preceding section.
%
% add utils folder to path (contains plotting functions)
addpath('utils')
% provide data directory and RUN
datadir = 'data';
%
%
%% (run and advance to next section)
%
channels_to_analyse = [4,5];
%
% DATA -------------------------------------------------------------------%
%
% DATASET INFO
%
% dataset name
load([datadir,slsh,'datasets.mat'],'datasets')
% channel label: YFP and RFP
load([datadir,slsh,'chlbls.mat'],'chlbls')
% weeks of chase
load([datadir,slsh,'time_points.mat'],'time_points')
% number os section in each dataset: confetti sample have only 1 section,
% kras samples have 3.
load([datadir,slsh,'store_sections.mat'],'store_sections')
%
% CLONAL AND CLUSTER DATA ------------------------------------------------%
%
% lobe areas
load([datadir,slsh,'lobe_area_total.mat'],'lobe_area_total')
% clone sizes (total no. of cells)
load([datadir,slsh,'clone_sizes_total.mat'],'clone_sizes_total')
% number of SPC+ cells in each clone
load([datadir,slsh,'clone_sizes_SPC_pos.mat'],'clone_sizes_SPC_pos')
%
% DATA FOR SIZE-DISTANCE ANALYSIS
%
% distance between pairs of RFP and YFP clusters in each dataset
load([datadir,slsh,'all_ds.mat'],'all_ds')
% total number of cells in the reference (central) RFP cluster
load([datadir,slsh,'all_nv_cent.mat'],'all_nv_cent')
% total number of cells in the surrounding (neighbouring) YFP cluster
load([datadir,slsh,'all_nv_neigh.mat'],'all_nv_neigh')
% number of spc+ cells in the surrounding (neighbouring) YFP cluster
load([datadir,slsh,'all_nv_neigh_spc.mat'],'all_nv_neigh_spc')
%
%-------------------------------------------------------------------------%
%
time_points_lbl = {}; for i=1:length(time_points); time_points_lbl{end+1} = num2str(time_points(i)); end

% MERGE DATA (run and advance to next section)
clone_sizes_combined = {}; clone_sizes_SPC_combined = {};
lobe_area_sum = {}; avg_clone_size = {}; lbls = {};
counter = 0;
for ndatas = 1:length(datasets)
    dataset = datasets{ndatas};
    for nmice = 1:length(clone_sizes_total{ndatas})
        lobe_area_sum{ndatas}(nmice,1) = sum(lobe_area_total{ndatas}{nmice}(:));
        if contains(dataset,'kras')
            % combine same mouse and same channel
            for nchannel = 1:length(channels_to_analyse)
                clone_sizes_combined{ndatas}{nmice,nchannel} = vertcat(clone_sizes_total{ndatas}{nmice}{:,nchannel});
                clone_sizes_SPC_combined{ndatas}{nmice,nchannel} = vertcat(clone_sizes_SPC_pos{ndatas}{nmice}{:,nchannel}); 
            end
        else
            % combine same mouse and both channels
            clone_sizes_combined{ndatas}{nmice,1} = vertcat(clone_sizes_total{ndatas}{nmice}{:});
            clone_sizes_SPC_combined{ndatas}{nmice,1} = vertcat(clone_sizes_SPC_pos{ndatas}{nmice}{:});
        end
    end

end

clone_number_density_per_lobe = {}; clone_number_per_lobe = {};
avg_clone_size = {}; spc_all = {}; lbls = {}; time_points_merge = [];
counter = 0;
for ndatas = 1:length(datasets)
    for nchannel = 1:size(clone_sizes_combined{ndatas},2)
        if size(clone_sizes_combined{ndatas},2) == 1
            lbls{end+1} = datasets{ndatas};
            time_points_merge(end+1,1) = time_points(ndatas);
            counter = counter + 1;
            for nmice = 1:length(clone_sizes_total{ndatas})
                cs = clone_sizes_combined{ndatas}{nmice,nchannel};
                nlobes = length(lobe_area_total{ndatas}{nmice});
                sections = store_sections{ndatas}(1);
                clone_number_per_lobe{counter}(nmice,:) = [length(cs)/nlobes/sections,sum(cs>1)/nlobes/sections,sum(cs==1)/nlobes/sections];
                clone_number_density_per_lobe{counter}(nmice,:) = [length(cs)/lobe_area_sum{ndatas}(nmice)/nlobes/sections,sum(cs>1)/lobe_area_sum{ndatas}(nmice)/nlobes/sections];
                avg_clone_size{counter}(nmice,:) = [mean(cs),mean(cs(cs>1))];
                spc = clone_sizes_SPC_combined{ndatas}{nmice,nchannel};
                spc_all{counter}(nmice,:) = [(sum(cs) - sum(spc))/lobe_area_sum{ndatas}(nmice)/nlobes/sections,sum(spc(cs>1))/lobe_area_sum{ndatas}(nmice)/nlobes/sections,sum(spc(cs>1))/nlobes/sections,sum(spc)/sum(cs)];

            end
        end
    end
end
for nchannel = 1:2
    for ndatas = 1:length(datasets)
        if size(clone_sizes_combined{ndatas},2) == 2
            lbls{end+1} = [datasets{ndatas},' ',chlbls{nchannel}];
            time_points_merge(end+1,1) = time_points(ndatas);
            counter = counter + 1;
            for nmice = 1:length(clone_sizes_total{ndatas})
                cs = clone_sizes_combined{ndatas}{nmice,nchannel};
                nlobes = length(lobe_area_total{ndatas}{nmice});
                sections = store_sections{ndatas}(1);
                clone_number_per_lobe{counter}(nmice,:) = [length(cs)/nlobes/sections,sum(cs>1)/nlobes/sections,sum(cs==1)/nlobes/sections];
                clone_number_density_per_lobe{counter}(nmice,:) = [length(cs)/lobe_area_sum{ndatas}(nmice)/nlobes/sections,sum(cs>1)/lobe_area_sum{ndatas}(nmice)/nlobes/sections];
                avg_clone_size{counter}(nmice,:) = [mean(cs),mean(cs(cs>1))];
                spc = clone_sizes_SPC_combined{ndatas}{nmice,nchannel};
                % | spc cell density | spc cell density (n>1)| |
                spc_all{counter}(nmice,:) = [(sum(cs) - sum(spc))/lobe_area_sum{ndatas}(nmice)/nlobes/sections,sum(spc(cs>1))/lobe_area_sum{ndatas}(nmice)/nlobes/sections,sum(spc(cs>1))/nlobes/sections,1-sum(spc)/sum(cs)];        
            end
        end
    end
end

time_points_mrg_lbl = {}; for i=1:length(time_points_merge); time_points_mrg_lbl{end+1} = num2str(time_points_merge(i)); end
time_points_merge(14:17) = time_points_merge(14:17)+5; 

%% FIGURE 1E: Combined CCDFs CONFETTI 

lbl_parcial = [];
cdfs = {};
ignore_clones_of_size = 1;

f = figure; f.Position = [000 100 350 320];
for ndatas = 1:length(datasets) %
    dataset = datasets{ndatas};
    for nchannel=1:size(clone_sizes_combined{ndatas},2)
        if ndatas<10
            for nmice= 1:size(clone_sizes_combined{ndatas},1)
                clone_sizes_tmp = clone_sizes_combined{ndatas}{nmice,nchannel};
                bincentsall = ignore_clones_of_size+0.5:1:10^5+0.5;
                if ~isempty(clone_sizes_tmp)
                    [cdf,bincents] = empirical_cdf(clone_sizes_tmp,0,nchannel,ignore_clones_of_size);
                    if length(cdf)<length(bincentsall); cdf(length(bincentsall)) = 0; end
                    if length(bincents)<length(bincentsall); bincents(length(bincentsall)) = 0; end
                    cdfs{ndatas,nchannel}(nmice,:) = cdf;
                end
            end
            c = cdfs{ndatas,nchannel};
            hold on
            cmean = mean(c);
            bincents_tmp = bincentsall;
            bincents_tmp(cmean==0) = [];
            cmean(cmean==0) = [];
            plot(bincents_tmp(1:end-1),cmean(1:end-1),'linewidth',1.5)
            lbl_parcial{end+1} = dataset;
        end
    end
end
box on;
xlabel('Clone size (cells)'); ylabel('Cumulative probability')
ylim([1e-4 1])
set(gca,'YScale','log'); set(gcf,'Color','w')
legend(lbl_parcial); legend boxoff

%% FIGURE S1D: LOBE AREA PER LOBE PER CONDITION
time_points_lbl_tmp = time_points_lbl;
time_points_lbl_tmp{2} = ''; time_points_lbl_tmp{3} = '';

f = figure; f.Position = [000 100 250 250];
    dot_plot_samples(lobe_area_total(1:9),time_points(1:9),time_points_lbl_tmp(1:9),'Lobe area (\mum^2)',-1); 
    xlim([(time_points(1))-5,(time_points(9))+5]); xlabel('Week')

%% % FIGURE S1K, S3P: Number of clones per lobe
time_points_lbl_tmp = time_points_lbl;
time_points_lbl_tmp{2} = '';

% FIGURE S2D
figure; 
dot_plot(clone_number_per_lobe(1:9),time_points_merge(1:9),time_points_lbl_tmp(1:9),'Clones count (>1 cell)',2);%,p_value_columns)
title('confetti')
xlim([(time_points_merge(1))-5,(time_points_merge(9))+5]); xlabel('Week'); yl = ylim(); ylim([0,yl(2)])

% FIGURE S3P
figure; 
dot_plot(clone_number_per_lobe(15:end),time_points_merge(15:end),time_points_mrg_lbl(15:end),'Clones count  (>1 cell)',2);%,p_value_columns)
title('kras')
xlim([(time_points_merge(15))-1,(time_points_merge(end))+1]); xlabel('Week'); yl = ylim(); ylim([0,yl(2)])

%% FIGURE S1L: Clone density CONFETTI
time_points_lbl_tmp = time_points_lbl;
time_points_lbl_tmp{2} = '';

f = figure; set(f,'position',[0,100,350,350])
dot_plot(clone_number_density_per_lobe(1:9),time_points_merge(1:9),time_points_mrg_lbl(1:9),'Clone density (count/\mum^2)',2);%,p_value_columns)
title('Confetti (n>1 cell)')
xlim([(time_points_merge(1))-3,(time_points_merge(9))+3]); xlabel('Week')
title('Confetti')


%% FIGURE S3Q: Clone density RFP+ KRAS
f = figure; set(f,'position',[0,100,600,250])
dot_plot(clone_number_density_per_lobe(15:end),time_points_merge(15:end),time_points_mrg_lbl(15:end),'Clone density (count/\mum^2)',2);%,p_value_columns)
xlim([(time_points_merge(15))-1,(time_points_merge(end))+1]); xlabel('Week')
title('RFP+ Kras')
%% FIGURE 5B: COMPARE CONFETTI AND WT AVG CLONE SIZE
f = figure; f.Position = [000 100 260 250];
[f,g,meanvals,xvals] = dot_plot(avg_clone_size(1:3),time_points_merge(1:3),time_points_mrg_lbl(1:3),'Average clone size (n)',2); 
fit_linear_growth(xvals,meanvals); 
hold on
t1 = 10;
[~,~,meanvals,xvals] = dot_plot(avg_clone_size(t1:13),time_points_merge(t1:13),time_points_mrg_lbl(t1:13),'Average clone size (>1 cell)',2);
[~,p1] = ttest2(avg_clone_size{1}(:,2),avg_clone_size{11}(:,2));
plot_p_value([time_points_merge(1)-0.1,time_points_merge(1)+0.1],5.8,0.05,round(p1,2));
[~,p1] = ttest2(avg_clone_size{2}(:,2),avg_clone_size{12}(:,2));
plot_p_value([time_points_merge(2)-0.1,time_points_merge(2)+0.1],6.4,0.05,round(p1,2));
[~,p1] = ttest2(avg_clone_size{3}(:,2),avg_clone_size{13}(:,2));
plot_p_value([time_points_merge(3)-0.1,time_points_merge(3)+0.1],6.4,0.05,round(p1,2));
hold off
% 
fit_linear_growth(xvals(1:2),meanvals(1:2)); 
fit_linear_growth(xvals(2:end),meanvals(2:end)); 
xlabel('Week') 
ylabel('Avg. Clone size (no. of cells)')
legend off
ylim([0,7])

%% FIGURE 1J, 2H, 2I, 2J, 6A, 6B, 6C, S2E, S2F, S3L 
% SPC+ vs total

counter = 0;
cs_all = {};
spc_all = {}; 
lbls = {};
for ndatas = 1:length(datasets)
    for nchannel = 1:size(clone_sizes_combined{ndatas},2)
        if size(clone_sizes_combined{ndatas},2) == 1
            lbls{end+1} = datasets{ndatas};
            counter = counter + 1;
            cs = []; spc = [];
            for nmice = 1:length(clone_sizes_total{ndatas})
                cs_tmp = clone_sizes_combined{ndatas}{nmice,nchannel};
                spc_tmp = clone_sizes_SPC_combined{ndatas}{nmice,nchannel}; 
                
                cs = [cs;cs_tmp];
                spc = [spc;spc_tmp];
            end
            cs_all{counter} = cs;
            spc_all{counter} = spc;
        end
    end
end
for nchannel = 1:2
    for ndatas = 1:length(datasets)
        if size(clone_sizes_combined{ndatas},2) == 2
            lbls{end+1} = [datasets{ndatas},' ',chlbls{nchannel}];
            counter = counter + 1;
            cs = []; spc = [];
            for nmice = 1:length(clone_sizes_total{ndatas})
                cs_tmp = clone_sizes_combined{ndatas}{nmice,nchannel};
                spc_tmp = clone_sizes_SPC_combined{ndatas}{nmice,nchannel}; 
                
                cs = [cs;cs_tmp];
                spc = [spc;spc_tmp];
            end
            cs_all{counter} = cs;
            spc_all{counter} = spc;
        end
    end
end
%
rm_clones_of_size = 1;
for i = 1:length(cs_all)
    f = figure; f.Position = [000 100 400 350];
    cs = cs_all{i};
    spc = spc_all{i};
    spc(cs <= rm_clones_of_size) = [];
    cs(cs <= rm_clones_of_size) = [];
    dsc = densityScatterChart(cs,100.*(1-spc./cs));
    pause(0.1)
    [tcl, ax, scat] = unmanage(dsc);

    cd = scat.CData; cd = 300*cd./max(cd);
        cd(cd<5) = 5;
    cd(cd>100) = 100;
    scat.SizeData = cd;  

    ylabel('Fraction of SPC- cells (%)')
    xlabel('Clone size (cells)')

    pause(1)
    title(lbls{i})

    a=colorbar;
    ylabel(a,'Clones')
    xl = xlim();
    if i >= 11 && i ~= 14
        set(ax,'xScale','log');
        xlim([1,10^ceil(log10(xl(2)))])
    else
        xlim([0,xl(2)])
    end
    ylim([0,100])
end

%% FIGURE S2B, S2C, S2D and MENDELEY Data Figure 1C, 1D, 1E, 1F, 1G, 1H 
% ratio of pro-Sftp^- cells / pro-Stfp^+ cells in RFP+ clones.

rm_clones_of_size = 1;

for i = 1:length(cs_all)
    f = figure; f.Position = [000 100 430 350];
    cs = cs_all{i};
    spc = spc_all{i};
    spc(cs <= rm_clones_of_size) = [];
    cs(cs <= rm_clones_of_size) = [];
    dsc = densityScatterChart(spc,cs-spc);
    pause(0.1)
    [tcl, ax, scat] = unmanage(dsc);

    ylabel('Number of SPC- cells')
    xlabel('Number of SPC+ cells')
    pause(1)
    title(lbls{i})
    a=colorbar;
    ylabel(a,'Clones')
    xl = xlim();
    yl = ylim();
    maxl = max([xl(2),yl(2)]); maxl = ceil(maxl/5); maxl = maxl*5;
    axis([0 maxl 0 maxl]);
end

%% FIGURE S1I (both panels)
cols = {[0.9290, 0.6940, 0.1250],[0.8500, 0.3250, 0.0980]};
mkrs = {'o','s'};

ignore_clones_of_size = 1;
for ndatas = [1,9] %
    f = figure;
    dataset = datasets{ndatas};
    for nchannel=1:2
        for nmice= 1:length(clone_sizes_total{ndatas})
            clone_sizes_tmp = vertcat(clone_sizes_total{ndatas}{nmice}{:,nchannel});
            hold on
            [emp_cdf,bincents] = empirical_cdf(clone_sizes_tmp,0,nchannel,ignore_clones_of_size);
            plot(bincents,emp_cdf,['-',mkrs{nchannel}],'color',cols{nchannel},'MarkerFaceColor','w')
            hold off
            set(gca,'YScale','log')
            set(gcf,'color','w');
            box on;

        end
        %             yl = ylim(); ylim([10e-4 1])
        %             xl = xlim(); xlim([0 xl(2)])
    end
    sgtitle(dataset)
end

%% ANALYSE DISRIBUTIONS OF CLONE SIZE
% HERE WE FIT THE BIEXPONENTIAL MODEL AND EXTRACT THE RELVANT PARAMETERS,
% f_F (here fs), ns (here nbar_p), nf (here nbar).
% RUN THIS SECTION BEFORE CONTINUING.

ratios = [];
fraction_of_SPC = {}; fs_all = {};
sigma_all = {}; Delta_s_all = {};
V_all = {}; nbar_p_all = {};
nbar_all = {};
all_cdfs = {};

ignore_clones_of_size = 1;
average_clone_sizes = {};
average_clone_sizes_nosing = {};
fit_model = 1;

show_plots = 0; % change to 1 to show plots with fits
counter = 0;
size_thesholds = [];

ratio = {};

for ndatas = 1:length(datasets) %

    dataset = datasets{ndatas};
    residuals = {};
    if show_plots
        f = figure;
        if size(clone_sizes_combined{ndatas},2) == 1
            f.Position = [000 100 300 300];
        else
            f.Position = [000 100 670 300];
        end
    end
    for nchannel=1:size(clone_sizes_combined{ndatas},2)
        counter = counter + 1;
        residuals{nchannel} = [];
        if ndatas>3 %&& ndatas~=10
            for nmice= 1:size(clone_sizes_combined{ndatas},1)

                clone_sizes_tmp = clone_sizes_combined{ndatas}{nmice,nchannel};
                clone_sizes_SPC_tmp = clone_sizes_SPC_combined{ndatas}{nmice,nchannel};
                clone_sizes_negSPC_tmp = clone_sizes_tmp - clone_sizes_SPC_tmp;

                if ~isempty(clone_sizes_tmp)
                    [cdf,bincents] = empirical_cdf(clone_sizes_tmp,0,nchannel,ignore_clones_of_size);
                    all_cdfs{ndatas}{nchannel}(nmice,:) = [cdf,zeros(1,100000-length(cdf))];
                    if numel(cdf) >= 5
                        if show_plots
                            subplot(1,size(clone_sizes_combined{ndatas},2),nchannel)
                            hold on

                            [emp_cdf,bincents] = empirical_cdf(clone_sizes_tmp,1,nchannel,ignore_clones_of_size);
                        else
                            [emp_cdf,bincents] = empirical_cdf(clone_sizes_tmp,0,nchannel,ignore_clones_of_size);
                        end
                        if max(clone_sizes_tmp)>1000
                            factor = 1/20;
                        else
                            factor = 0.5;
                        end

                        % we ignore the largest clone in the fit
                        if ndatas >= 11; clone_sizes_tmp(clone_sizes_tmp == max(clone_sizes_tmp))=[]; end

                        size_threshold_min = minimize(clone_sizes_tmp,ignore_clones_of_size,factor);
                        size_thesholds{ndatas}(nmice,nchannel) = size_threshold_min;
                        if size_threshold_min>1

                            if ndatas == 11 && nmice == 3 % in this sample there is a break in the distribution. We prioritize a good fit of smaller clones
                                upper_cuttoff = 70;
                                [gof,nbar,nbar_p,n0int,T1p,~,zero_all,zero_long] = fit_biexponential(clone_sizes_tmp,size_threshold_min,ignore_clones_of_size,upper_cuttoff);
                            else
                                [gof,nbar,nbar_p,n0int,T1p,~,zero_all,zero_long] = fit_biexponential(clone_sizes_tmp,size_threshold_min,ignore_clones_of_size);
                            end
                            ratio{ndatas}(nmice,nchannel) = zero_all/zero_long;

                            Cn = theoretical_cdf(nbar,nbar_p,n0int,ignore_clones_of_size,0,max(bincents));

                            if show_plots
                                hold on
                                plot(bincents,Cn,'--k','HandleVisibility','off')
                                hold off

                                xlabel('Clone size (cells)'); ylabel('Cumulative probability')
                            end
                            r = 1;
                            [Delta_s,sigma,fs] = estimate_remaining_parameters(r,time_points(ndatas),nbar,n0int,T1p);
                            nbar_all{counter}(nmice,1) = nbar;
                            nbar_p_all{counter}(nmice,1) = nbar_p;
                            Delta_s_all{counter}(nmice,nchannel) = Delta_s;
                            sigma_all{counter}(nmice,1) = sigma;
                            fs_all{counter}(nmice,1) = fs;
                            residuals{nchannel}{end+1} = [bincents',emp_cdf',Cn'];
                        else
                            residuals{nchannel}{end+1} = [];
                        end

                        if show_plots
                            set(gca,'YScale','log')
                            if contains(dataset,'kras')
                                title([chlbls{nchannel},' Clones'])
                            else
                                title('Combined Clones')
                            end
                            set(gcf,'color','w');
                            box on;
                        end
                    end
                end
            end
            if show_plots
                yl = ylim(); ylim([10e-4 1])
                xl = xlim(); xlim([0 xl(2)])
            end
        end
    end
end

%% FIGURE 1G, 1H, 1I, 2F, 2G, S6D, S6E, S6F
% Previous section must by executed before generating these plots

figure;
% CONFETTI
subplot(3,3,3); 
[~,~,meanvals,xvals] = dot_plot(fs_all(4:9),time_points_merge(4:9),time_points_mrg_lbl(4:9),'f_F',1);
h = refline(0,mean(meanvals)); h.Color = [0,0,0]; h.LineStyle = ':';
xlabel('Week post induction')
xlim([(time_points_merge(4))-3,(time_points_merge(9))+3]);
legend({'data mean\pmSD','linear fit'}); legend boxoff;
% CONFETTI
subplot(3,3,2);
title('Confetti')
[~,~,meanvals,xvals] = dot_plot(nbar_all(4:9),time_points_merge(4:9),time_points_mrg_lbl(4:9),'nbar',1);
fit_linear_growth(time_points_merge(4:9),meanvals); 
xlabel('Week post induction'); ylabel('$\overline{n}_f$','Interpreter','Latex')
xlim([(time_points_merge(4))-3,(time_points_merge(9))+3]);
yl = ylim();
ylim([0 ceil(yl(2))]);
legend({'data mean\pmSD','linear fit'},'Location','northwest'); legend boxoff;
% CONFETTI
subplot(3,3,1); 
[~,~,meanvals,xvals] = dot_plot(nbar_p_all(4:9),time_points_merge(4:9),time_points_mrg_lbl(4:9),'nbar_p',1);
fit_linear_growth(time_points_merge(4:9),meanvals);
yl = ylim();
ylim([0 yl(2)]);
xlabel('Week post induction'); ylabel('$\overline{n}_s$','Interpreter','Latex')
xlim([(time_points_merge(4))-3,(time_points_merge(9))+3]);
yl = ylim();
ylim([0 ceil(yl(2))]);

indata = [10,12,14,16];
intts = 10:13;
legend({'data mean\pmSD','linear fit'},'Location','northwest'); legend boxoff;
% sgtitle('Confetti')

% KRAS (WT)
subplot(3,3,6); 
[~,~,meanvals,xvals] = dot_plot(fs_all(indata),time_points_merge(intts),time_points_mrg_lbl(intts),'f_F',1);
h = refline(0,mean(meanvals)); h.Color = [0,0,0]; h.LineStyle = ':';
xlabel('Week post induction')
legend({'data','average'}); legend boxoff;
% KRAS (WT)
subplot(3,3,5); 
[~,~,meanvals,xvals,variance_weight] = dot_plot(nbar_all(indata),time_points_merge(intts),time_points_mrg_lbl(intts),'nbar',1);
fit_linear_growth(time_points_merge(intts(1:2)),meanvals(1:2));
fit_linear_growth(time_points_merge(intts(2:end)),meanvals(2:end));
xlabel('Week post induction'); ylabel('$\overline{n}_f$','Interpreter','Latex')
title('Kras (WT)')
yl = ylim();
ylim([0 30]);
legend({'data mean\pmSD','linear fit'},'Location','southeast'); legend boxoff;
% KRAS (WT)
subplot(3,3,4); 
[~,~,meanvals,xvals,variance_weight] = dot_plot(nbar_p_all(indata),time_points_merge(intts),time_points_mrg_lbl(intts),'nbar_p',1);
fit_linear_growth(time_points_merge(intts(1:2)),meanvals(1:2));
fit_linear_growth(time_points_merge(intts(2:end)),meanvals(2:end));
xlabel('Week post induction'); ylabel('$\overline{n}_s$','Interpreter','Latex')
yl = ylim();
ylim([0 5]);
indata = [11,13,15,17];
intts = 14:17;
legend({'data mean\pmSD','linear fit'},'Location','southeast'); legend boxoff;

% KRAS (RFP+)
subplot(3,3,9); 
[~,~,meanvals,xvals] = dot_plot(fs_all(indata),time_points_merge(intts)-time_points_merge(min(intts))+1,time_points_mrg_lbl(intts),'f_F',1);
h = refline(0,mean(meanvals)); h.Color = [0,0,0]; h.LineStyle = ':';
xlabel('Week post induction')
ylim([0 0.15]);
disp('Warning: ylim fixed manually')
legend({'data mean\pmSD','linear fit'}); legend boxoff;

% KRAS (RFP+)
subplot(3,3,8); 
[~,~,meanvals,xvals] = dot_plot(nbar_all(indata),time_points_merge(intts)-time_points_merge(min(intts))+1,time_points_mrg_lbl(intts),'nbar',1);
fit_exponential_growth(xvals(end-1:end),meanvals(end-1:end));
fit_exponential_growth(xvals(1:end-1),meanvals(1:end-1));
title('Kras (RFP+)') 
xlabel('Week post induction'); ylabel('$\overline{n}_f$','Interpreter','Latex')
legend({'data mean\pmSD','exp fit'},'Location','northwest'); legend boxoff;
set(gca,'yscale','log')

% KRAS (RFP+)
subplot(3,3,7); 
[~,~,meanvals,xvals] = dot_plot(nbar_p_all(indata),time_points_merge(intts)-time_points_merge(min(intts))+1,time_points_mrg_lbl(intts),'nbar_p',1);
fit_exponential_growth(xvals(end-1:end),meanvals(end-1:end));
fit_exponential_growth(xvals(1:end-1),meanvals(1:end-1));
xlabel('Week post induction'); ylabel('$\overline{n}_s$','Interpreter','Latex')
yl = ylim();
ylim([0 yl(2)]);
set(gca,'yscale','log')
set(gcf,'position',[10,10,700,800])
legend({'data mean\pmSD','exp fit'},'Location','northwest'); legend boxoff;

%% Related to FIGURE 1F, 2E, 5D, 5E, S1N, S1O, S3I, S3J, S6C
% PLOT CCDFs OF CLONE SIZE (average +- SD)
% This shows the experimental curves only, for plots with model run 
% simulation script "sim_two_pop_model.m"

for ndatas = 4:length(datasets) %
    for chan = 1:size(clone_sizes_combined{ndatas},2)
        f = figure;
        dataset = datasets{ndatas};
        cdfs = all_cdfs{ndatas}{chan};

        rec_cdfvar = var(cdfs,1);
        cdf = mean(cdfs,1);
        x = ignore_clones_of_size+1:length(cdf)+1;

        rec_cdfvar(cdf == 0) = [];
        x(cdf == 0) = [];
        cdf(cdf == 0) = [];
        errorbar(x,cdf,sqrt(rec_cdfvar),'.')
        set(gca,'yscale','log')
        set(gcf,'Color','w')
        ylim([0.99e-5,1])
        f.Position = [100 200 350 350];
        if size(clone_sizes_combined{ndatas},2) == 1
            title(dataset)
        else
            title([dataset,' ',chlbls{chan}])
        end
        xlabel('Clone size (cells)'); ylabel('Cumulative probability')
    end
end

%% FIGURE 5G, 5H, 5I, 5J, 6E, 6F, S6G, S6H, S6I
% Size-size correlation

channels_to_analyse = {}; chcombo = {};
channels_to_analyse{end+1} = [5,4]; chcombo{end+1} = 'RFP-YFP'; chlbls = {'RFP','YFP'};
distsbin = 0:50:1500;
for ndatas = 1:length(datasets)
    dataset = datasets{ndatas};
    for nchannel = 1:length(channels_to_analyse)

        nv_cent = all_nv_cent{ndatas,nchannel};
        nv_neigh = all_nv_neigh{ndatas,nchannel};
        nv_neigh_spc = all_nv_neigh_spc{ndatas,nchannel};
        ds = all_ds{ndatas,nchannel};

        nv_neigh_spc = nv_neigh_spc./nv_neigh;
        [~,~,bins] = histcounts(ds,distsbin);

        bincent = distsbin(2:end);
        pos = bincent(bins);
        und = unique(pos);
        f = figure; f.Position = [000 100 800 250];

        subplot(1,3,1)
        alld = {}; for i = 1:length(pos); alld{end+1} = ['<',num2str(pos(i))]; end
        ord = {};
        avgs = [];

        avg_cent = [];
        avg_spcneg_frac_neigh = []; std_spcneg_frac_neigh = [];
        for i = 1:length(und)
            ord{end+1} = ['<',num2str(und(i))];
            avgs(i) = mean(nv_neigh(pos == und(i)));
            avg_cent(i) = mean(nv_cent(pos == und(i)));
            nv_neigh_spc_neg = 1 - nv_neigh_spc;
            avg_spcneg_frac_neigh(i) = mean(nv_neigh_spc_neg(pos == und(i)));
            std_spcneg_frac_neigh(i) = std(nv_neigh_spc_neg(pos == und(i)));
        end
        vs = violinplot(nv_neigh, alld,'GroupOrder',ord,'ShowMean',true);
        hold on
        plot(1:length(und),avgs,'k')
        hold off
        yl = ylim();
        ylim([0,yl(2)])

        set(gcf,'Color','w')
        xlabel('Dist. to Kras+ clone (um)')
        ylabel('WT Clone size (n)')

        subplot(1,3,2)
        sublbl = {};
        hold on
        for i = 1:1:6
            in = find(pos == und(i)); n = nv_neigh(in); binedges = 0.5:1:max(n)+0.5; bincents = (binedges(2:end) + binedges(1:end - 1))/2; Nnorm = histcounts(n,binedges,'Normalization','pdf');
            cdf = 1 - cumsum(Nnorm)./sum(Nnorm);
            bincents(cdf<10e-6) = [];
            cdf(cdf<10e-6) = [];

            sublbl{end+1} = ['d < ',num2str(und(i)),' \mum'];
            plot(bincents,cdf,'LineWidth',1)
            set(gca,'yScale','log')
        end
        hold off
        box on
        xlabel('WT Clone size (cells)')
        ylabel('Cumulative probability')
        legend(sublbl); legend boxoff
        set(gcf,'Color','w')
        yl = ylim();
        ylim([0,yl(2)])

        subplot(1,3,3)
        vs = violinplot(100*(1-nv_neigh_spc), alld,'GroupOrder',ord,'ShowMean',true);
        set(gcf,'Color','w')
        xlabel('Dist. to Kras+ clone (um)')
        ylabel('SPC- cells in WT clones (%)')
        yl = ylim();
        ylim([0,yl(2)])
        xl = xlim();
        xlim([0,xl(2)])

        hold on
        errorbar(1:length(und),100*avg_spcneg_frac_neigh,std_spcneg_frac_neigh,'k')
        hold off
        yl = ylim();
        ylim([0,100])
        xl = xlim();
        xlim([0,xl(2)])
        sgtitle(dataset)
    end
end

%% FUNCTIONS

function [f,g,meanvals,xvals,variance_weight] = dot_plot(data,x_ticks,x_labels,y_label,clumns,sem) 
meanvals = [];
variance_weight = [];
xvals = [];

hold on
f = []; g = [];
for i = 1:length(data)
    if clumns>0
        dat = data{i}(:,clumns);
    else
        dat = data{i}(:);
    end

    dat = vertcat(dat(:));
    if nargin == 6
        norm = sqrt(numel(dat));
    else
        norm = 1;
    end
    if numel(dat) == 1
        f(i) = scatter(x_ticks(i),dat,'filled','MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'HandleVisibility','off');
        g(i) = errorbar(x_ticks(i),dat,0,'sk','LineWidth',1,'MarkerFaceColor','k','HandleVisibility','off');
        meanvals(i) = dat;
        variance_weight(i) = 1;
    else
        f(i) = scatter(x_ticks(i)*ones(size(dat))+(rand(size(dat))-0.5)/5,dat,'filled','MarkerEdgeColor','k','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2,'HandleVisibility','off');
        if i == length(data)
            g(i) = errorbar(x_ticks(i),mean(dat),std(dat)/norm,'sk','LineWidth',1,'MarkerFaceColor','k');
        else
            g(i) = errorbar(x_ticks(i),mean(dat),std(dat)/norm,'sk','LineWidth',1,'MarkerFaceColor','k','HandleVisibility','off');
        end
        meanvals(i) = mean(dat);
        variance_weight(i) = var(dat);
        xvals(i) = x_ticks(i);
    end
end
hold off
xticks(x_ticks)
xticklabels(x_labels)
ylabel(y_label)
xlim([0,x_ticks(end)+1])
box on
set(gcf,'Color','w')
end

function [meanvals,xvals] = dot_plot_samples(data,x_ticks,x_labels,y_label,clumns) % data is a cell array
meanvals = [];
xvals = [];
hold on
for i = 1:length(data)
    dat_all = data{i};
    for nmouse = 1:length(dat_all)
        dat = dat_all{nmouse};
        [c,m] = candm(nmouse);
        if clumns>0
            dat = dat(:,clumns);
        end

        if numel(dat) == 1
            scatter(x_ticks(i),dat,m,'MarkerEdgeColor',c);
        else
            scatter(x_ticks(i).*ones(size(dat))+(rand(size(dat))-0.5)/2,dat,'filled',m,'MarkerEdgeColor',c);
        end
    end
    d = vertcat(dat_all{:});
    if clumns>0
     d = d(:,clumns);
    end
    errorbar(x_ticks(i),mean(d),std(d),'sk','LineWidth',1,'MarkerFaceColor','k')
    meanvals(i) = mean(d);
    xvals(i) = x_ticks(i);
end
hold off
xticks(x_ticks)
xticklabels(x_labels)
ylabel(y_label)
box on
set(gcf,'Color','w')
end

function plot_p_value(x,y,yfactor,p)
hold on
if isstring(p)
    text(mean(x),y,p,'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    plot(x,[y,y],'-k','HandleVisibility','off')
%     %    ylim([0.5,length(data)+0.5])
%     axis([0.5,length(data)+0.5 0 maxy*1.2])
else
    text(mean(x),y+yfactor,['p = ',num2str(p)],'HorizontalAlignment','Center','VerticalAlignment','Bottom')
    plot(x,[y,y],'-k','HandleVisibility','off')
end
hold off
end

function [c,m] = candm(n)
% define colours and markers for plots
cols = {[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980],[0.9290, 0.6940, 0.1250],[0.4940, 0.1840, 0.5560],[0.4660, 0.6740, 0.1880],[0.3010, 0.7450, 0.9330],[0.6350, 0.0780, 0.1840]};
mkrs = {'o','s','d','*','^','v'};
c = cols{n};
m = mkrs{n};
end

function [size_threshold_min, gof] = minimize(clone_sizes,ignore_clones_of_size,factor)
% Input:
% clone_sizes: list of clone sizes for a given sample
% find optimal size_threshold_min in the range [1,30] cells. This range was
% chosen as it was obseved than in all samples the transition between the
% two exponential regimes lied within this range.
opts = optimset('FunValCheck','on','MaxFunEvals',000,'MaxIter',3000,'TolX',1e-5);
lb = 3;
ub = max(clone_sizes)*factor;
if max(clone_sizes)*factor<=lb; ub=max(clone_sizes); end
[size_threshold_min, gof] = fminbnd(@(size_threshold) fit_biexponential(clone_sizes,size_threshold,ignore_clones_of_size),lb,ub,opts);
end

function [gof,nbar,nbar_p,p0p,T1,bincents,zero_all,zero_long] = fit_biexponential(clone_sizes,size_threshold,ignore_clones_of_size,upper_cuttoff)
% ---------------------------------------------------------------------
% This function fits the bi-exponential CDF function to the empirical
% data as described in the "Model fits" section of the METHODS S1 PDF file
% ---------------------------------------------------------------------
% construct empirical CDF
[emp_cdf,bincents] = empirical_cdf(clone_sizes,0,1,ignore_clones_of_size);
% extract values above size_threshold (long term dependence)
if nargin < 4
    upper_cuttoff = inf;
end

cdf_long = emp_cdf(bincents > size_threshold & bincents < upper_cuttoff);
bincents_long = bincents(bincents > size_threshold & bincents < upper_cuttoff);
% ignore zeros in the CDF
bincents_long(cdf_long == 0) = []; cdf_long(cdf_long == 0) = [];

% define exponential function for fitting
exp_fun = @(c, xi) c(1)*exp(-xi./c(2));

% 1) fit exponential decay to the long term dependence
% initial guess
x0 = [1,mean(clone_sizes)];
% run fitting routine
% provide reasonable lower (lb) and upper (ub) bounds to improve
% convergence
lb = [0 1]; ub = [10 max(clone_sizes)];
c_long = lsqcurvefit(exp_fun, x0, bincents_long, cdf_long,lb,ub);
% evaluate exponential function with the fitted parameters
y_long = exp_fun(c_long,bincents_long);
% ---------------------------------------------------------------------
% extract nbar
nbar = c_long(2);
% extract the composite parameter p0p=fs*(2r-1)/r) from the n=0 intersect
p0p = c_long(1);
% ---------------------------------------------------------------------
% 2) Substract the long-term dependence from the empirical CDF and
% estimate the short term decay
y_long = exp_fun(c_long,bincents);
cdf_short = emp_cdf - y_long;
% ignore zeros in the CDF

cdf_short = cdf_short(bincents <= size_threshold);
bincents_short = bincents(bincents <= size_threshold);
bincents_short = bincents_short(cdf_short > 0); cdf_short(cdf_short <= 0) = [];
if ~isempty(bincents_short)
    % extract values below size_threshold (short term dependence)
    % initial guess
    x0 = [1,size_threshold];
    % run fitting routine
    % provide reasonable lower (lb) and upper (ub) bounds to improve
    % convergence
    lb = [0.5 1]; ub = [1 max(clone_sizes)];
    [c_short] = lsqcurvefit(exp_fun, x0, bincents_short, cdf_short,lb,ub);
    % evaluate exponential function with the fitted parameters

    y_short = exp_fun(c_short,bincents_short);
    T1 = 1-y_short(1); % y_short(1) corresponds to the Prob. that a clone has > than 1 cell.

    % ---------------------------------------------------------------------
    % extract nbar_p
    nbar_p = c_short(2);
    % ---------------------------------------------------------------------
    % evaluate the theoretical cdf, and extract the goodness-of-fit
    % parameter S.
    Cn = theoretical_cdf(nbar,nbar_p,p0p,ignore_clones_of_size,0,max(clone_sizes));
    emp_cdf(emp_cdf<0) = 0;
    [R2,S] = gof_estimation(bincents(1:end-1),log(emp_cdf(1:end-1)),log(Cn(1:end-1)));
    gof = S;
    [emp_cdfz,bc] = empirical_cdf(clone_sizes,0,1,0);
    zero_all = emp_cdfz(1);
    zero_long = exp_fun(c_long,1);
else
    gof = 10^10;
end
end

function [R2,S] = gof_estimation(bincents,emp_cdf,theo_cdf)
    % we just compare clone sizes where we do have experimental data
    theo_cdf = theo_cdf(isfinite(emp_cdf));
    emp_cdf = emp_cdf(isfinite(emp_cdf));
    ybar = nanmean(emp_cdf);
    Stot = nansum((emp_cdf - ybar).^2);
    Sres = nansum((emp_cdf - theo_cdf).^2);
    R2 = 1 - Sres/Stot;
    S = sqrt(Sres/length(bincents));
end

function Cn = theoretical_cdf(nbar,nbar_p,n0int,ignore_clones_of_size,show_plot,max_cs)
    % Here we contuct the theoretical cdf based on the input parameters
    % nbar, nbar_p and p0p.
    % If ignore_clones_of_size > 0, then T{n<=ignore_clones_of_size} is removed from Tn(t).
    % provide range of clone sizes to evaluate
    bincents = 0:1:10^5;
    Tn = @(xi) n0int*exp(-xi./nbar)./nbar+(1-n0int)*exp(-xi./nbar_p)./nbar_p;
    Tn_list = Tn(bincents)/(1-sum(Tn(0:ignore_clones_of_size)));
    Cn = 1 - cumsum(Tn_list)./sum(Tn_list);
    Cn = Cn(ignore_clones_of_size+1:max_cs);
    if show_plot == 1
       plot(bincents(ignore_clones_of_size+1:max_cs),Cn,'--k','HandleVisibility','off') 
    end
end

function [cdf,bincents] = empirical_cdf(clone_sizes,show_plots,nchannel,ignore_clones_of_size)
% we contuct the empirical cdf
% define bins (1 cell spacing)
binedges = (ignore_clones_of_size+0.5):1:max(clone_sizes)*1+0.5; bincents = (binedges(2:end) + binedges(1:end - 1))/2;
% compute pdf
Nnorm = histcounts(clone_sizes,binedges,'Normalization','pdf');
% compute cdf
cdf = 1 - cumsum(Nnorm)./sum(Nnorm);
if show_plots == 1
    [c,m] = candm(nchannel);
    plot(bincents(1:end-1),cdf(1:end-1),m,'color',c,'MarkerFaceColor','w','MarkerSize',5)
end
end % empirical_cdf

function [fitresult, gof] = fit_exponential_growth(x,y)
[xData, yData] = prepareCurveData( x, y );
ft = fittype( 'exp1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.00525470190820502 1.47524019659865];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
hold on
plot( fitresult ,':k');
end

function [fitresult, gof] = fit_linear_growth(x,y,er)
if nargin < 3
    er = [ones(1,length(x))];
end
[xData, yData, weights] = prepareCurveData( x, y, er );
% Set up fittype and options.
ft = fittype( 'poly1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Weights = weights;
opts.StartPoint = [0.00525470190820502 1.47524019659865];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
hold on
plot( fitresult ,':k');
end

function [Delta_s,sigma,fs] = estimate_remaining_parameters(r,t0,nbar,n0int,T1p)
% calculate parameters for a given value of "r":
factor = (2*r-1)/r;
% stem cell expansion rate
Delta_s = (1/t0)*log(factor^2*nbar);
% stem cell cycling rate
sigma = Delta_s/(2*r-1);
% T1p being the n=1 intersect of the small clone exponential decay
fs = n0int*(1-T1p)/factor; % from fsp = fs/(1-s(t))
end

% DATASETS
function [dir_path,file_name,nsections,nlobe] = data_to_analyse(dataset,root_path)

if ispc; slsh = '\'; else; slsh = '/'; end

dir_path = {};
file_name = {};
nsections = [];
nlobe = [];
% CORN OIL CONTROL
switch dataset
    case {'cornoil','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'corn oil controls',slsh,'JLAL22.3c lobe1 merged tile scan z stack'];
        file_name{end+1} = 'JLAL22.3c lobe1 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [0 1 1 1]; % [condition | timepoint | sample | lobe]
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'corn oil controls',slsh,'JLAL22.3c lobe2 merged tile scan z stack'];
        file_name{end+1} = 'JLAL22.3c lobe2 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [0 1 1 2];
end
switch dataset
    case {'conf1w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4a merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4a merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 1 1]; % [condition | timepoint | sample | lobe]
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4a merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4a merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4a merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.4a merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4a merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.4a merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 1 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4a merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL22.4a merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 1 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4b merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4b merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4b merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4b merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4b merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.4b merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4b merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.4b merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 2 4];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4g merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4g merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4g merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4g merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4g merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.4g merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'1week',slsh,'new dose',slsh,'JLAL22.4g merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.4g merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 1 3 4];
end
switch dataset
    case {'conf2w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2a lobe1 merged tile scan z stack'];
        file_name{end+1} = 'JLAL30.2a lobe1 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2a lobe2 merged tile scan z stack'];
        file_name{end+1} = 'JLAL30.2a lobe2 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2a lobe3 merged tile scan z stack'];
        file_name{end+1} = 'JLAL30.2a lobe3 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2a lobe4 merged tile scan z stack'];
        file_name{end+1} = 'JLAL30.2a lobe4 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 1 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2a lobe5 merged tile scan z stack'];
        file_name{end+1} = 'JLAL30.2a lobe5 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 1 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2b merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL30.2b merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2b merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL30.2b merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2b merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL30.2b merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2b merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL30.2b merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 2 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL30.2b merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL30.2b merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 2 5];

        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL36.2e lobe1 merged tile scan z stack'];
        file_name{end+1} = 'JLAL36.2e lobe1 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL36.2e lobe2 merged tile scan z stack'];
        file_name{end+1} = 'JLAL36.2e lobe2 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL36.2e lobe3 merged tile scan z stack'];
        file_name{end+1} = 'JLAL36.2e lobe3 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL36.2e lobe4 merged tile scan z stack'];
        file_name{end+1} = 'JLAL36.2e lobe4 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'2weeks',slsh,'JLAL36.2e lobe5 merged tile scan z stack'];
        file_name{end+1} = 'JLAL36.2e lobe5 merged tile scan z stack.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 2 3 5];
end
switch dataset
    case {'conf4w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4c merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4c merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 4 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4c merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4c merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 4 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4c merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.4c merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 4 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4c merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.4c merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 4 1 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4c merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL22.4c merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1;  nlobe(end+1,:) = [1 4 1 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4d merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4d merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4d merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4d merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 2 2];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4h merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4h merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4h merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4h merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4h merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.4h merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4h merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.4h merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL22.4h merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL22.4h merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 3 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL29.1b merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL29.1b merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 4 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL29.1b merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL29.1b merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 4 2];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL29.1d merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL29.1d merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 5 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'4weeks',slsh,'new dose',slsh,'JLAL29.1d merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL29.1d merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 4 5 2];
end
switch dataset
    case {'conf8w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'8weeks',slsh,'low dose',slsh,'JLAL24.2j merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.2j merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 8 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'8weeks',slsh,'low dose',slsh,'JLAL24.2j merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.2j merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 8 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'8weeks',slsh,'low dose',slsh,'JLAL24.2j merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.2j merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 8 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'8weeks',slsh,'low dose',slsh,'JLAL24.2j merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL24.2j merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 8 1 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'8weeks',slsh,'low dose',slsh,'JLAL24.2j merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL24.2j merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 8 1 5];
        %
end
switch dataset
    case {'conf12w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1a merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.1a merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1a merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.1a merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1a merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.1a merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1a merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL24.1a merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 1 4];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1e merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.1e merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1e merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.1e merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1e merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.1e merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1e merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL24.1e merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 2 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL24.1e merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL24.1e merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 2 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL25.1a merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL25.1a merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL25.1a merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL25.1a merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL25.1a merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL25.1a merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL25.1a merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL25.1a merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'12weeks',slsh,'new dose',slsh,'JLAL25.1a merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL25.1a merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 12 3 5];
end
switch dataset
    case {'conf24w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1c merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.1c merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1c merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.1c merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1c merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.1c merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 1 3];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1g merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.1g merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1g merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.1g merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1g merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.1g merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1g merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL24.1g merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 2 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL24.1g merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL24.1g merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 2 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL25.1c merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL25.1c merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL25.1c merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL25.1c merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL25.1c merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL25.1c merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL25.1c merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL25.1c merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'24weeks',slsh,'low dose',slsh,'JLAL25.1c merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL25.1c merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 24 3 5];
end
switch dataset
    case {'conf36w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5a merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL23.5a merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5a merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL23.5a merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5a merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL23.5a merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5a merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL23.5a merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 1 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5a merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL23.5a merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 1 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5g merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL23.5g merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5g merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL23.5g merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5g merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL23.5g merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL23.5g merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL23.5g merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 2 4];
        %
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL24.3k lobe1 merged tile scan z stack'];
        file_name{end+1} = 'JLAL24.3k lobe1 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL24.3k lobe2 merged tile scan z stack'];
        file_name{end+1} = 'JLAL24.3k lobe2 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL24.3k lobe3 merged tile scan z stack'];
        file_name{end+1} = 'JLAL24.3k lobe3 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL24.3k lobe4 merged tile scan z stack'];
        file_name{end+1} = 'JLAL24.3k lobe4 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'36weeks',slsh,'JLAL24.3k lobe5 merged tile scan z stack'];
        file_name{end+1} = 'JLAL24.3k lobe5 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 36 3 5];
end
switch dataset
    case {'conf52w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL22.4e merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.4e merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL22.4e merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.4e merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL22.4e merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.4e merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL22.4e merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.4e merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 1 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL22.4e merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL22.4e merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 1 5];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5h merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL23.5h merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5h merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL23.5h merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5h merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL23.5h merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5h merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL23.5h merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 2 4];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5i merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL23.5i merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5i lobe2 merged tile scan z stack'];
        file_name{end+1} = 'JLAL23.5i lobe2 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5i merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL23.5i merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5i lobe4 merged tile scan z stack'];
        file_name{end+1} = 'JLAL23.5i lobe4 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'52weeks',slsh,'JLAL23.5i lobe5 merged tile scan z stack'];
        file_name{end+1} = 'JLAL23.5i lobe5 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 52 3 5];
        %

end
switch dataset
    case {'conf60w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6c merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.6c merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6c merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.6c merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6c merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.6c merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6c merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.6c merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 1 4];
        %
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6d merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.6d merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6d merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.6d merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6d merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.6d merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6d merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.6d merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 2 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6d merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL22.6d merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 2 5];

        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6f merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL22.6f merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6f merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL22.6f merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6f merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL22.6f merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6f merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL22.6f merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'60weeks',slsh,'JLAL22.6f merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL22.6f merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 60 3 5];
        %
end
switch dataset
    case {'conf72w','allconf','all'}
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2b merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.2b merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 1 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2b merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.2b merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 1 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2b merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.2b merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 1 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2b merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL24.2b merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 1 4];
        
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2h merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL24.2h merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 2 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2h merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL24.2h merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 2 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2h merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL24.2h merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 2 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL24.2h merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL24.2h merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 2 4];
        %
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1c merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL28.1c merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 3 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1c merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL28.1c merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 3 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1c merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL28.1c merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 3 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1c merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL28.1c merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 3 4];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1c merged tile scan z stack lobe5'];
        file_name{end+1} = 'JLAL28.1c merged tile scan z stack lobe5.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 3 5];

        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1d merged tile scan z stack lobe1'];
        file_name{end+1} = 'JLAL28.1d merged tile scan z stack lobe1.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 4 1];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1d merged tile scan z stack lobe2'];
        file_name{end+1} = 'JLAL28.1d merged tile scan z stack lobe2.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 4 2];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1d merged tile scan z stack lobe3'];
        file_name{end+1} = 'JLAL28.1d merged tile scan z stack lobe3.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 4 3];
        dir_path{end+1} = [root_path,slsh,'confetti',slsh,'72weeks',slsh,'JLAL28.1d merged tile scan z stack lobe4'];
        file_name{end+1} = 'JLAL28.1d merged tile scan z stack lobe4.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [1 72 4 4];
        %
end

switch dataset
    case {'kras4d','allkras','all'}
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG47.3a lobe1 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG47.3a lobe1 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 1 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG47.3a lobe2 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG47.3a lobe2 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 1 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG47.3a lobe3 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG47.3a lobe3 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 1 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG47.3a lobe4 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG47.3a lobe4 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 1 4];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG47.3a lobe5 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG47.3a lobe5 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 1 5];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG78.1b lobe1 merged tile scan z stacks.lif'];
        file_name{end+1} = 'JLAG78.1b lobe1 merged tile scan z stacks.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 2 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG78.1b lobe2 merged tile scan z stacks.lif'];
        file_name{end+1} = 'JLAG78.1b lobe2 merged tile scan z stacks.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 2 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG78.1b lobe3 merged tile scan z stacks.lif'];
        file_name{end+1} = 'JLAG78.1b lobe3 merged tile scan z stacks.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 2 3];
        
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG78.1h lobe3 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG78.1h lobe3 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 3 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4days',slsh,'JLAG78.1h lobe4 merged tile scan z stack.lif'];
        file_name{end+1} = 'JLAG78.1h lobe4 merged tile scan z stack.lif.tif'; nsections(end+1) = 1; nlobe(end+1,:) = [2 4*0.1429 3 2];
end

switch dataset
    case {'kras1w','allkras','all'}
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1a lobe1 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1a lobe1 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 1 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1a lobe2 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1a lobe2 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 1 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1a lobe3 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1a lobe3 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 1 3];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1b lobe1 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1b lobe1 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 2 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1b lobe3 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1b lobe3 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 2 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1b lobe4 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1b lobe4 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 2 3];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1d merged lobe1 tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1d merged lobe1 tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 3 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1d merged lobe2 tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1d merged lobe2 tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 3 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1d lobe3 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1d lobe3 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 3 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1d merged lobe6 tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1d merged lobe6 tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 3 4];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1e lobe1 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1e lobe1 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 4 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1e lobe2 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1e lobe2 merged tile scan z stacks.lif - stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 4 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1e lobe5 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1e lobe5 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 4 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'1week',slsh,'JLAG38.1e lobe6 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG38.1e lobe6 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 1 4 4];
end

switch dataset
    case {'kras2w','allkras','all'}
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1a lobe2 merged tile scans'];
        file_name{end+1} = 'JLAG30.1a lobe2 merged tile scans.lif - Stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 1 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1a lobe3 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1a lobe3 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 1 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1a lobe4 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1a lobe4 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 1 3];
        
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1b lobe2 merged tile scan z stack'];
        file_name{end+1} = 'JLAG30.1b lobe2 merged tile scan z stack.lif - Stack '; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 2 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1b lobe5 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1b lobe5 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 2 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1b lobe6 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1b lobe6 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 2 4];
                
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1d lobe1 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1d lobe1 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 3 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1d lobe5 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1d lobe5 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 3 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG30.1d lobe6 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG30.1d lobe6 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 3 3];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG68.3a lobe2 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG68.3a lobe2 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 4 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG68.3a lobe4 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG68.3a lobe4 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 4 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG68.3a lobe6 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG68.3a lobe6 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 4 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG68.3a lobe7 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG68.3a lobe7 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 4 4];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'2weeks',slsh,'JLAG68.3a lobe8 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG68.3a lobe8 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 2 4 5];
        % %
end

switch dataset
    case {'kras4w','allkras','all'}
%         % % %% kras 4 weeks
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG31.1a lobe1 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG31.1a lobe1 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 1 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG31.1a lobe5 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG31.1a lobe5 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 1 2];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG31.1a lobe6 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG31.1a lobe6 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 1 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG31.1a lobe8 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG31.1a lobe8 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 1 5];
      
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3a lobe2 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3a lobe2 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 2 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3a lobe3 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3a lobe3 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 2 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3a lobe7 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3a lobe7 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 2 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3a lobe8 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3a lobe8 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 2 4];

        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3g lobe1 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3g lobe1 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 3 1];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3g lobe2 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3g lobe2 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 3 2];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3g lobe3 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3g lobe3 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 3 3];
        dir_path{end+1} = [root_path,slsh,'kras',slsh,'4weeks',slsh,'JLAG53.3g lobe4 merged tile scan z stacks'];
        file_name{end+1} = 'JLAG53.3g lobe4 merged tile scan z stacks.lif - stack'; nsections(end+1) = 3; nlobe(end+1,:) = [2 4 3 4];
end

end

