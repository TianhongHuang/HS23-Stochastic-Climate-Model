%--------------------------------------------------------------------------
% [SBR]: HISTOGRAM
% [TAG]: histogram data
% [INTRO]: comparing of histogram data of rainfall statistics in
%          standard simulation v.s climate-change study
%          sensitivity studies of cloud albedo A_b/A_f
%--------------------------------------------------------------------------
% @author:  Tianhong Huang
% @date:    May 2023   
% @program: SH23
% @version: v5.1
%--------------------------------------------------------------------------

%% LOADING OPS FROM SIMULATIONS
%--------------------------------------------------------------------------
% please choose the right captions in Ab or Af test
%--------------------------------------------------------------------------
str_caption_1 = 'A_b -0.1';
str_caption_2 = 'A_b -0.05';
str_caption_3 = 'standard';
str_caption_4 = 'A_b +0.05';
str_caption_5 = 'A_b +0.1';

%str_caption_1 = 'A_f -0.1';
%str_caption_2 = 'A_f -0.05';
%str_caption_3 = 'standard';
%str_caption_4 = 'A_f +0.05';
%str_caption_5 = 'A_f +0.1';

str_caption_11 = 'standard';
str_caption_12 = 'Increased CO_2';

%--------------------------------------------------------------------------
% 'hist_rain_events.mat' includes histogram data of rainfall events size
% + rainfall events duration + dry spell duration
%--------------------------------------------------------------------------
% 'hist_cluster_size.mat' includes histogram data of cluster size of both
% shallow stratocumulus cloud and deep convective cloud
%--------------------------------------------------------------------------
% 'hist_xxx_1.mat' = simulation with albedo value -0.1
%--------------------------------------------------------------------------
load('hist_rain_events_1.mat')
rshist_1 = rshist; % rain-size-histogram data
rthist_1 = rthist; % rain-time-histogram data
dthist_1 = dthist; % dry-time-histogram data
sthist_1 = sthist; % shallowcluster-time-histogram data
%--------------------------------------------------------------------------
% xxx_re = normalized data (rescaled)
%--------------------------------------------------------------------------
rshist_1_re = rshist_re; 
rthist_1_re = rthist_re;
dthist_1_re = dthist_re;
sthist_1_re = sthist_re;
load('hist_cluster_size_1.mat')
scshist_1 = scs_hist; % shallow-cluster-size-histogram data
dcshist_1 = dcs_hist; % deep-convection-size-histogram data
scshist_1_re = scs_hist_rescale;
dcshist_1_re = dcs_hist_rescale;
%--------------------------------------------------------------------------
% 'hist_xxx_2.mat' = simulation with albedo value -0.05
%--------------------------------------------------------------------------
load('hist_rain_events_2.mat')
rshist_2 = rshist;
rthist_2 = rthist;
dthist_2 = dthist;
sthist_2 = sthist;
rshist_2_re = rshist_re;
rthist_2_re = rthist_re;
dthist_2_re = dthist_re;
sthist_2_re = sthist_re;
load('hist_cluster_size_2.mat')
scshist_2 = scs_hist;
dcshist_2 = dcs_hist;
scshist_2_re = scs_hist_rescale;
dcshist_2_re = dcs_hist_rescale;
%--------------------------------------------------------------------------
% 'hist_xxx_3.mat' = simulation with albedo value +0.0
%--------------------------------------------------------------------------
load('hist_rain_events_3.mat')
rshist_3 = rshist;
rthist_3 = rthist;
dthist_3 = dthist;
sthist_3 = sthist;
rshist_3_re = rshist_re;
rthist_3_re = rthist_re;
dthist_3_re = dthist_re;
sthist_3_re = sthist_re;
load('hist_cluster_size_3.mat')
scshist_3 = scs_hist;
dcshist_3 = dcs_hist;
scshist_3_re = scs_hist_rescale;
dcshist_3_re = dcs_hist_rescale;
%--------------------------------------------------------------------------
% 'hist_xxx_4.mat' = simulation with albedo value +0.05
%--------------------------------------------------------------------------
load('hist_rain_events_4.mat')
rshist_4 = rshist;
rthist_4 = rthist;
dthist_4 = dthist;
sthist_4 = sthist;
rshist_4_re = rshist_re;
rthist_4_re = rthist_re;
dthist_4_re = dthist_re;
sthist_4_re = sthist_re;
load('hist_cluster_size_4.mat')
scshist_4 = scs_hist;
dcshist_4 = dcs_hist;
scshist_4_re = scs_hist_rescale;
dcshist_4_re = dcs_hist_rescale;
%--------------------------------------------------------------------------
% 'hist_xxx_5.mat' = simulation with albedo value +0.1
%--------------------------------------------------------------------------
load('hist_rain_events_5.mat')
rshist_5 = rshist;
rthist_5 = rthist;
dthist_5 = dthist;
sthist_5 = sthist;
rshist_5_re = rshist_re;
rthist_5_re = rthist_re;
dthist_5_re = dthist_re;
sthist_5_re = sthist_re;
load('hist_cluster_size_5.mat')
scshist_5 = scs_hist;
dcshist_5 = dcs_hist;
scshist_5_re = scs_hist_rescale;
dcshist_5_re = dcs_hist_rescale;
%--------------------------------------------------------------------------
% 'hist_xxx_11.mat' = standard simulation
%--------------------------------------------------------------------------
load('hist_rain_events_11.mat')
rshist_11 = rshist;
rthist_11 = rthist;
dthist_11 = dthist;
sthist_11 = sthist;
rshist_11_re = rshist_re;
rthist_11_re = rthist_re;
dthist_11_re = dthist_re;
sthist_11_re = sthist_re;
load('hist_cluster_size_11.mat')
scshist_11 = scs_hist;
dcshist_11 = dcs_hist;
scshist_11_re = scs_hist_rescale;
dcshist_11_re = dcs_hist_rescale;
%--------------------------------------------------------------------------
% 'hist_xxx_12.mat' = climate-change simulation
%--------------------------------------------------------------------------
load('hist_rain_events_12.mat')
rshist_12 = rshist;
rthist_12 = rthist;
dthist_12 = dthist;
sthist_12 = sthist;
rshist_12_re = rshist_re;
rthist_12_re = rthist_re;
dthist_12_re = dthist_re;
sthist_12_re = sthist_re;
load('hist_cluster_size_12.mat')
scshist_12 = scs_hist;
dcshist_12 = dcs_hist;
scshist_12_re = scs_hist_rescale;
dcshist_12_re = dcs_hist_rescale;


%% PLOTTING
%--------------------------------------------------------------------------
% save plot figures in folder [output_fig]
%--------------------------------------------------------------------------
cd([pwd,'/output_fig/'])

%--------------------------------------------------------------------------
% FIG_1: histogram of rain events (sens-study)
%--------------------------------------------------------------------------
fig_rainhist_ss = figure;
set(gcf,'units','points','position',[10,10,1200,400])

subplot(1,2,1)
line11 = loglog(rsbinscale, rshist_1_re/sum(rshist_1_re),'--d','linewidth',1.5);
hold on
line12 = loglog(rsbinscale, rshist_2_re/sum(rshist_2_re),'--o','linewidth',1.5);
line13 = loglog(rsbinscale, rshist_3_re/sum(rshist_3_re),'-s','linewidth',1.5);
line14 = loglog(rsbinscale, rshist_4_re/sum(rshist_4_re),'-.+','linewidth',1.5);
line15 = loglog(rsbinscale, rshist_5_re/sum(rshist_5_re),'-.*','linewidth',1.5);
set(line11, 'MarkerFaceColor', get(line11,'Color'));
set(line12, 'MarkerFaceColor', get(line12,'Color'));
set(line13, 'MarkerFaceColor', get(line12,'Color'));
set(line14, 'MarkerFaceColor', get(line12,'Color'));
set(line15, 'MarkerFaceColor', get(line12,'Color'));
title('A. rain events size','FontWeight','Bold','FontSize',18)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim('auto')
xlabel('rainfall size (mm)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)
line21 = loglog(rtbinscale, rthist_1_re/sum(rthist_1_re),'--d','linewidth',1.5);
hold on
line22 = loglog(rtbinscale, rthist_2_re/sum(rthist_2_re),'--o','linewidth',1.5);
line23 = loglog(rtbinscale, rthist_3_re/sum(rthist_3_re),'-s','linewidth',1.5);
line24 = loglog(rtbinscale, rthist_4_re/sum(rthist_4_re),'-.+','linewidth',1.5);
line25 = loglog(rtbinscale, rthist_5_re/sum(rthist_5_re),'-.*','linewidth',1.5);
set(line21, 'MarkerFaceColor', get(line11,'Color'));
set(line22, 'MarkerFaceColor', get(line21,'Color'));
set(line23, 'MarkerFaceColor', get(line23,'Color'));
set(line24, 'MarkerFaceColor', get(line24,'Color'));
set(line25, 'MarkerFaceColor', get(line25,'Color'));
title('B. rain events duration','FontWeight','Bold','FontSize',18)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim('auto')
xlabel('rainfall time (mins)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_rainhist_ss,'[sens][hist]rainfall');
close(fig_rainhist_ss)

%--------------------------------------------------------------------------
% FIG_2: histogram of convective cluster size (sens-study)
%--------------------------------------------------------------------------
fig_dschist_ss = figure;

line31 = loglog(csbinscale(2:end)*25,dcshist_1_re(2:end)/sum(dcshist_1_re),'--d','linewidth',1.5);
hold on
line32 = loglog(csbinscale(2:end)*25,dcshist_2_re(2:end)/sum(dcshist_2_re),'--o','linewidth',1.5);
line33 = loglog(csbinscale(2:end)*25,dcshist_3_re(2:end)/sum(dcshist_3_re),'-s','linewidth',1.5);
line34 = loglog(csbinscale(2:end)*25,dcshist_4_re(2:end)/sum(dcshist_4_re),'-.+','linewidth',1.5);
line35 = loglog(csbinscale(2:end)*25,dcshist_5_re(2:end)/sum(dcshist_5_re),'-.*','linewidth',1.5);
set(line31, 'MarkerFaceColor', get(line31,'Color'));
set(line32, 'MarkerFaceColor', get(line32,'Color'));
set(line33, 'MarkerFaceColor', get(line33,'Color'));
set(line34, 'MarkerFaceColor', get(line34,'Color'));
set(line35, 'MarkerFaceColor', get(line35,'Color'));
title('deep cluster size histogram','FontWeight','Bold','FontSize',18)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim('auto')
xlabel('cluster size (km^2)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_dschist_ss,'[sens][hist]conv_cs');
close(fig_dschist_ss)

%--------------------------------------------------------------------------
% FIG_3: histogram of dry events (sens_study)
%--------------------------------------------------------------------------
fig_dthist_ss = figure;

line41 = loglog(dtbinscale, dthist_1_re,'--d','linewidth',1.5);
hold on
line42 = loglog(dtbinscale, dthist_2_re,'--o','linewidth',1.5);
line43 = loglog(dtbinscale, dthist_3_re,'-s','linewidth',1.5);
line44 = loglog(dtbinscale, dthist_4_re,'-.+','linewidth',1.5);
line45 = loglog(dtbinscale, dthist_5_re,'-.*','linewidth',1.5);
set(line41, 'MarkerFaceColor', get(line41,'Color'));
set(line42, 'MarkerFaceColor', get(line41,'Color'));
set(line43, 'MarkerFaceColor', get(line43,'Color'));
set(line44, 'MarkerFaceColor', get(line44,'Color'));
set(line45, 'MarkerFaceColor', get(line45,'Color'));
title('dry events duration histogram','FontWeight','Bold','FontSize',18)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim('auto')
xlabel('dry time (mins)','FontWeight','Bold','FontSize',16)
ylabel('number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_dthist_ss,'[sens][hist]drytime');
close(fig_dthist_ss)

%--------------------------------------------------------------------------
% FIG_4: histogram of shallow cloud (sens-study)
%--------------------------------------------------------------------------
fig_schist_ss = figure;

set(gcf,'units','points','position',[10,10,1200,400])

subplot(1,2,1)
line51 = loglog(stbinscale, sthist_1_re/sum(sthist_1_re),'--d','linewidth',1.5);
hold on
line52 = loglog(stbinscale, sthist_2_re/sum(sthist_2_re),'--o','linewidth',1.5);
line53 = loglog(stbinscale, sthist_3_re/sum(sthist_3_re),'-s','linewidth',1.5);
line54 = loglog(stbinscale, sthist_4_re/sum(sthist_4_re),'-.+','linewidth',1.5);
line55 = loglog(stbinscale, sthist_5_re/sum(sthist_5_re),'-.*','linewidth',1.5);
set(line51, 'MarkerFaceColor', get(line51,'Color'));
set(line52, 'MarkerFaceColor', get(line52,'Color'));
set(line53, 'MarkerFaceColor', get(line53,'Color'));
set(line54, 'MarkerFaceColor', get(line54,'Color'));
set(line55, 'MarkerFaceColor', get(line55,'Color'));
title('A. shallow cluster duration','FontWeight','Bold','FontSize',18)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim('auto')
xlabel('cluster duration (mins)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)
line61 = loglog(csbinscale(2:end)*25,scshist_1_re(2:end)/sum(scshist_1_re),'--d','linewidth',1.5);
hold on
line62 = loglog(csbinscale(2:end)*25,scshist_2_re(2:end)/sum(scshist_2_re),'--o','linewidth',1.5);
line63 = loglog(csbinscale(2:end)*25,scshist_3_re(2:end)/sum(scshist_3_re),'-s','linewidth',1.5);
line64 = loglog(csbinscale(2:end)*25,scshist_4_re(2:end)/sum(scshist_4_re),'-.+','linewidth',1.5);
line65 = loglog(csbinscale(2:end)*25,scshist_5_re(2:end)/sum(scshist_5_re),'-.*','linewidth',1.5);
set(line61, 'MarkerFaceColor', get(line61,'Color'));
set(line62, 'MarkerFaceColor', get(line62,'Color'));
set(line63, 'MarkerFaceColor', get(line63,'Color'));
set(line64, 'MarkerFaceColor', get(line64,'Color'));
set(line65, 'MarkerFaceColor', get(line65,'Color'));
title('B. shallow cluster size','FontWeight','Bold','FontSize',18)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim('auto')
xlabel('cluster size (km^2)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_schist_ss,'[sens][hist]shallow_cluster');
close(fig_schist_ss)

%--------------------------------------------------------------------------
% FIG_5: histogram of rain events (climate-change)
%--------------------------------------------------------------------------
fig_rainhist_cc = figure;
set(gcf,'units','points','position',[10,10,1200,400])

subplot(1,2,1)
line11 = loglog(rsbinscale, rshist_11_re/sum(rshist_11_re),'-s','linewidth',1.5);
hold on
line12 = loglog(rsbinscale, rshist_12_re/sum(rshist_12_re),'--d','linewidth',1.5);
set(line11, 'MarkerFaceColor', get(line11,'Color'));
set(line12, 'MarkerFaceColor', get(line12,'Color'));
title('A. rain events size','FontWeight','Bold','FontSize',18)
legend(str_caption_11,str_caption_12,'location','best')
xlim('auto')
xlabel('rainfall size (mm)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)
line21 = loglog(rtbinscale, rthist_11_re/sum(rthist_11_re),'-s','linewidth',1.5);
hold on
line22 = loglog(rtbinscale, rthist_12_re/sum(rthist_12_re),'--d','linewidth',1.5);
set(line21, 'MarkerFaceColor', get(line21,'Color'));
set(line22, 'MarkerFaceColor', get(line22,'Color'));
title('B. rain events duration','FontWeight','Bold','FontSize',18)
legend(str_caption_11,str_caption_12,'location','best')
xlim('auto')
xlabel('rainfall time (mins)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_rainhist_cc,'[climate-change][hist]rainfall');
close(fig_rainhist_cc)

%--------------------------------------------------------------------------
% FIG_6: histogram of convective cluster size (climate-change)
%--------------------------------------------------------------------------
fig_dschist_cc = figure;

line31 = loglog(csbinscale(2:end)*25,dcshist_11_re(2:end)/sum(dcshist_11_re),'-s','linewidth',1.5);
hold on
line32 = loglog(csbinscale(2:end)*25,dcshist_12_re(2:end)/sum(dcshist_12_re),'--d','linewidth',1.5);
set(line31, 'MarkerFaceColor', get(line31,'Color'));
set(line32, 'MarkerFaceColor', get(line32,'Color'));
title('deep cluster size histogram','FontWeight','Bold','FontSize',18)
legend(str_caption_11,str_caption_12,'location','best')
xlim('auto')
xlabel('cluster size (km^2)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_dschist_cc,'[climate-change][hist]conv_cs');
close(fig_dschist_cc)

%--------------------------------------------------------------------------
% FIG_7: histogram of dry events (climate-change)
%--------------------------------------------------------------------------
fig_dthist_cc = figure;

line41 = loglog(dtbinscale, dthist_11_re,'-s','linewidth',1.5);
hold on
line42 = loglog(dtbinscale, dthist_12_re,'--d','linewidth',1.5);
set(line41, 'MarkerFaceColor', get(line41,'Color'));
set(line42, 'MarkerFaceColor', get(line42,'Color'));
title('dry events duration histogram','FontWeight','Bold','FontSize',18)
legend(str_caption_11,str_caption_12,'location','best')
xlim('auto')
xlabel('dry time (mins)','FontWeight','Bold','FontSize',16)
ylabel('number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_dthist_cc,'[climate-change][hist]drytime');
close(fig_dthist_cc)

%--------------------------------------------------------------------------
% FIG_8: histogram of shallow cloud (climate-change)
%--------------------------------------------------------------------------
fig_schist_cc = figure;

set(gcf,'units','points','position',[10,10,1200,400])

subplot(1,2,1)
line51 = loglog(stbinscale, sthist_11_re/sum(sthist_11_re),'-s','linewidth',1.5);
hold on
line52 = loglog(stbinscale, sthist_12_re/sum(sthist_12_re),'--d','linewidth',1.5);
set(line51, 'MarkerFaceColor', get(line51,'Color'));
set(line52, 'MarkerFaceColor', get(line52,'Color'));
title('A. shallow cluster duration','FontWeight','Bold','FontSize',18)
legend(str_caption_11,str_caption_12,'location','best')
xlim('auto')
xlabel('cluster duration (mins)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)
line61 = loglog(csbinscale(2:end)*25,scshist_11_re(2:end)/sum(scshist_11_re),'-s','linewidth',1.5);
hold on
line62 = loglog(csbinscale(2:end)*25,scshist_12_re(2:end)/sum(scshist_12_re),'--d','linewidth',1.5);
set(line61, 'MarkerFaceColor', get(line61,'Color'));
set(line62, 'MarkerFaceColor', get(line62,'Color'));
title('B. shallow cluster size','FontWeight','Bold','FontSize',18)
legend(str_caption_11,str_caption_12,'location','best')
xlim('auto')
xlabel('cluster size (km^2)','FontWeight','Bold','FontSize',16)
ylabel('normalized number of events','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_schist_cc,'[climate-change][hist]shallow_cluster');
close(fig_schist_cc)

cd ..
