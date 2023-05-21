%--------------------------------------------------------------------------
% [SBR]: SENS_STUDY
% [TAG]: sensitivity studies
% [INTRO]: comparing of time-meridionally-average and time series of
%          model variables in sensitivity studies of cloud albeos A_b/A_f
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

%--------------------------------------------------------------------------
% 'stats_varbs_1.mat' = simulation with albedo value -0.1
%--------------------------------------------------------------------------
load('stats_varbs_1.mat')
To_tavg_1 = To_tavg;
Tb_tavg_1 = Tb_tavg;
Tf_tavg_1 = Tf_tavg;
v1_tavg_1 = v1_tavg;
qf_tavg_1 = qf_tavg;
qb_tavg_1 = qb_tavg;
scloud_tavg_1 = scloud_tavg;
dcloud_tavg_1 = dcloud_tavg;
bcloud_tavg_1 = bcloud_tavg;
To_savg_1 = To_savg;
Tb_savg_1 = Tb_savg;
Tf_savg_1 = Tf_savg;
qf_savg_1 = qf_savg;
qb_savg_1 = qtb_savg;
scloud_savg_1 = scloud_rate_savg;
dcloud_savg_1 = dcloud_rate_savg;
bcloud_savg_1 = bcloud_rate_savg;
v1_yavg_1 = mean(v1_tavg_1);
To_yavg_1 = mean(To_tavg_1);
Tb_yavg_1 = mean(Tb_tavg_1);
Tf_yavg_1 = mean(Tf_tavg_1);
qf_yavg_1 = mean(qf_tavg_1);
qb_yavg_1 = mean(qb_tavg_1);
sigc_yavg_1 = mean(scloud_tavg_1);
sigf_yavg_1 = mean(dcloud_tavg_1);
sigb_yavg_1 = mean(bcloud_tavg_1);
%--------------------------------------------------------------------------
% 'stats_varbs_2.mat' = simulation with albedo value -0.05
%--------------------------------------------------------------------------
load('stats_varbs_2.mat')
To_tavg_2 = To_tavg;
Tb_tavg_2 = Tb_tavg;
Tf_tavg_2 = Tf_tavg;
v1_tavg_2 = v1_tavg;
qf_tavg_2 = qf_tavg;
qb_tavg_2 = qb_tavg;
scloud_tavg_2 = scloud_tavg;
dcloud_tavg_2 = dcloud_tavg;
bcloud_tavg_2 = bcloud_tavg;
To_savg_2 = To_savg;
Tb_savg_2 = Tb_savg;
Tf_savg_2 = Tf_savg;
qf_savg_2 = qf_savg;
qb_savg_2 = qtb_savg;
scloud_savg_2 = scloud_rate_savg;
dcloud_savg_2 = dcloud_rate_savg;
bcloud_savg_2 = bcloud_rate_savg;
v1_yavg_2 = mean(v1_tavg_2);
To_yavg_2 = mean(To_tavg_2);
Tb_yavg_2 = mean(Tb_tavg_2);
Tf_yavg_2 = mean(Tf_tavg_2);
qf_yavg_2 = mean(qf_tavg_2);
qb_yavg_2 = mean(qb_tavg_2);
sigc_yavg_2 = mean(scloud_tavg_2);
sigf_yavg_2 = mean(dcloud_tavg_2);
sigb_yavg_2 = mean(bcloud_tavg_2);
%--------------------------------------------------------------------------
% 'stats_varbs_3.mat' = simulation with albedo value +0.0
%--------------------------------------------------------------------------
load('stats_varbs_3.mat')
To_tavg_3 = To_tavg;
Tb_tavg_3 = Tb_tavg;
Tf_tavg_3 = Tf_tavg;
v1_tavg_3 = v1_tavg;
qf_tavg_3 = qf_tavg;
qb_tavg_3 = qb_tavg;
scloud_tavg_3 = scloud_tavg;
dcloud_tavg_3 = dcloud_tavg;
bcloud_tavg_3 = bcloud_tavg;
To_savg_3 = To_savg;
Tb_savg_3 = Tb_savg;
Tf_savg_3 = Tf_savg;
qf_savg_3 = qf_savg;
qb_savg_3 = qtb_savg;
scloud_savg_3 = scloud_rate_savg;
dcloud_savg_3 = dcloud_rate_savg;
bcloud_savg_3 = bcloud_rate_savg;
v1_yavg_3 = mean(v1_tavg_3);
To_yavg_3 = mean(To_tavg_3);
Tb_yavg_3 = mean(Tb_tavg_3);
Tf_yavg_3 = mean(Tf_tavg_3);
qf_yavg_3 = mean(qf_tavg_3);
qb_yavg_3 = mean(qb_tavg_3);
sigc_yavg_3 = mean(scloud_tavg_3);
sigf_yavg_3 = mean(dcloud_tavg_3);
sigb_yavg_3 = mean(bcloud_tavg_3);
%--------------------------------------------------------------------------
% 'stats_varbs_4.mat' = simulation with albedo value +0.05
%--------------------------------------------------------------------------
load('stats_varbs_4.mat')
To_tavg_4 = To_tavg;
Tb_tavg_4 = Tb_tavg;
Tf_tavg_4 = Tf_tavg;
v1_tavg_4 = v1_tavg;
qf_tavg_4 = qf_tavg;
qb_tavg_4 = qb_tavg;
scloud_tavg_4 = scloud_tavg;
dcloud_tavg_4 = dcloud_tavg;
bcloud_tavg_4 = bcloud_tavg;
To_savg_4 = To_savg;
Tb_savg_4 = Tb_savg;
Tf_savg_4 = Tf_savg;
qf_savg_4 = qf_savg;
qb_savg_4 = qtb_savg;
scloud_savg_4 = scloud_rate_savg;
dcloud_savg_4 = dcloud_rate_savg;
bcloud_savg_4 = bcloud_rate_savg;
v1_yavg_4 = mean(v1_tavg_4);
To_yavg_4 = mean(To_tavg_4);
Tb_yavg_4 = mean(Tb_tavg_4);
Tf_yavg_4 = mean(Tf_tavg_4);
qf_yavg_4 = mean(qf_tavg_4);
qb_yavg_4 = mean(qb_tavg_4);
sigc_yavg_4 = mean(scloud_tavg_4);
sigf_yavg_4 = mean(dcloud_tavg_4);
sigb_yavg_4 = mean(bcloud_tavg_4);
%--------------------------------------------------------------------------
% 'stats_varbs_5.mat' = simulation with albedo value +0.1
%--------------------------------------------------------------------------
load('stats_varbs_5.mat')
To_tavg_5 = To_tavg;
Tb_tavg_5 = Tb_tavg;
Tf_tavg_5 = Tf_tavg;
v1_tavg_5 = v1_tavg;
qf_tavg_5 = qf_tavg;
qb_tavg_5 = qb_tavg;
scloud_tavg_5 = scloud_tavg;
dcloud_tavg_5 = dcloud_tavg;
bcloud_tavg_5 = bcloud_tavg;
To_savg_5 = To_savg;
Tb_savg_5 = Tb_savg;
Tf_savg_5 = Tf_savg;
qf_savg_5 = qf_savg;
qb_savg_5 = qtb_savg;
scloud_savg_5 = scloud_rate_savg;
dcloud_savg_5 = dcloud_rate_savg;
bcloud_savg_5 = bcloud_rate_savg;
v1_yavg_5 = mean(v1_tavg_5);
To_yavg_5 = mean(To_tavg_5);
Tb_yavg_5 = mean(Tb_tavg_5);
Tf_yavg_5 = mean(Tf_tavg_5);
qf_yavg_5 = mean(qf_tavg_5);
qb_yavg_5 = mean(qb_tavg_5);
sigc_yavg_5 = mean(scloud_tavg_5);
sigf_yavg_5 = mean(dcloud_tavg_5);
sigb_yavg_5 = mean(bcloud_tavg_5);

xaxis = 5: 5: 10000;
timeaxis = data_time/365;

%% PLOTTING
%--------------------------------------------------------------------------
% save plot figures in folder [output_fig]
%--------------------------------------------------------------------------
cd([pwd,'/output_fig/'])

%--------------------------------------------------------------------------
% FIG_1: comparison of time-moridionally-avged model variables
%--------------------------------------------------------------------------
fig_yavg_varbs = figure;
set(gcf,'units','points','position',[10,10,1200,1600])

% meridional avg To
subplot(8,1,1)
plot(xaxis, To_yavg_1, '--', 'linewidth', 1.5)
hold on
plot(xaxis, To_yavg_2, '--', 'linewidth', 1.5)
plot(xaxis, To_yavg_3, '-', 'linewidth', 1.5)
plot(xaxis, To_yavg_4, '-.', 'linewidth', 1.5)
plot(xaxis, To_yavg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('T_o(K)','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('A. Ocean temperature','Fontweight','Bold','Fontsize',18)

% meridional avg Tb
subplot(8,1,2)
plot(xaxis, Tb_yavg_1, '--', 'linewidth', 1.5)
hold on
plot(xaxis, Tb_yavg_2, '--', 'linewidth', 1.5)
plot(xaxis, Tb_yavg_3, '-', 'linewidth', 1.5)
plot(xaxis, Tb_yavg_4, '-.', 'linewidth', 1.5)
plot(xaxis, Tb_yavg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('T_b(K)','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('B. Boundary layer temperature','Fontweight','Bold','Fontsize',18)

% meridional avg qb
subplot(8,1,3)
plot(xaxis, qb_yavg_1*1e+3, '--', 'linewidth', 1.5)
hold on
plot(xaxis, qb_yavg_2*1e+3, '--', 'linewidth', 1.5)
plot(xaxis, qb_yavg_3*1e+3, '-', 'linewidth', 1.5)
plot(xaxis, qb_yavg_4*1e+3, '-.', 'linewidth', 1.5)
plot(xaxis, qb_yavg_5*1e+3, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('q_b(mm)','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('C. Column water vapor, boundary layer','Fontweight','Bold','Fontsize',18)

% meridional avg shallow cloud
subplot(8,1,4)
plot(xaxis, sigb_yavg_1, '--', 'linewidth', 1.5)
hold on
plot(xaxis, sigb_yavg_2, '--', 'linewidth', 1.5)
plot(xaxis, sigb_yavg_3, '-', 'linewidth', 1.5)
plot(xaxis, sigb_yavg_4, '-.', 'linewidth', 1.5)
plot(xaxis, sigb_yavg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','west')
ylabel('\sigma_s','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
ylim([0.25 1])
title('D. Shallow cloud fraction','Fontweight','Bold','Fontsize',18)

% meridional avg Tf
subplot(8,1,5)
plot(xaxis, Tf_yavg_1, '--', 'linewidth', 1.5)
hold on
plot(xaxis, Tf_yavg_2, '--', 'linewidth', 1.5)
plot(xaxis, Tf_yavg_3, '-', 'linewidth', 1.5)
plot(xaxis, Tf_yavg_4, '-.', 'linewidth', 1.5)
plot(xaxis, Tf_yavg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('T_f(K)','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('E. Free troposphere temperature','Fontweight','Bold','Fontsize',18)

% meridional avg qf
subplot(8,1,6)
plot(xaxis, qf_yavg_1*1e+3, '--', 'linewidth', 1.5)
hold on
plot(xaxis, qf_yavg_2*1e+3, '--', 'linewidth', 1.5)
plot(xaxis, qf_yavg_3*1e+3, '-', 'linewidth', 1.5)
plot(xaxis, qf_yavg_4*1e+3, '-.', 'linewidth', 1.5)
plot(xaxis, qf_yavg_5*1e+3, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('q_f(mm)','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('F. Column water vapor, Free troposphere','Fontweight','Bold','Fontsize',18)

% meridional avg deep cloud
subplot(8,1,7)
plot(xaxis, sigf_yavg_1, '--', 'linewidth', 1.5)
hold on
plot(xaxis, sigf_yavg_2, '--', 'linewidth', 1.5)
plot(xaxis, sigf_yavg_3, '-', 'linewidth', 1.5)
plot(xaxis, sigf_yavg_4, '-.', 'linewidth', 1.5)
plot(xaxis, sigf_yavg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('\sigma_f','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('G. Deep cloud fraction','Fontweight','Bold','Fontsize',18)

% meridional avg zonal wind velocity
subplot(8,1,8)
plot(xaxis, v1_yavg_1, '--', 'linewidth', 1.5)
hold on
plot(xaxis, v1_yavg_2, '--', 'linewidth', 1.5)
plot(xaxis, v1_yavg_3, '-', 'linewidth', 1.5)
plot(xaxis, v1_yavg_4, '-.', 'linewidth', 1.5)
plot(xaxis, v1_yavg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,str_caption_3,str_caption_4,str_caption_5,'location','best')
ylabel('u_1 (m/s)','FontSize',16)
write_axis_x()
ax = gca;
ax.FontSize = 16;
title('H. First baroclinic mode of zonal wind velocity','Fontweight','Bold','Fontsize',18)

savefig(fig_yavg_varbs,'[sens][TMA]varbs.fig');
close(fig_yavg_varbs)

cd ..

%--------------------------------------------------------------------------
% FIG_2: comparison of time-series of model variables
%--------------------------------------------------------------------------
cd([pwd,'/output_fig/'])

fig_savg_varbs = figure;
set(gcf,'units','points','position',[10,10,1200,1000])

% To
subplot(5,1,1)
plot(timeaxis, To_savg_1, '--', 'linewidth', 1.5)
hold on
plot(timeaxis, To_savg_2, '--', 'linewidth', 1.5)
plot(timeaxis, To_savg_3, '-', 'linewidth', 1.5)
plot(timeaxis, To_savg_4, '-.', 'linewidth', 1.5)
plot(timeaxis, To_savg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1, str_caption_2, str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim([0,5])
ylabel('T_o(K)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('A. Ocean temperature, domain averaged','FontWeight','Bold','FontSize',18)

% Tb
subplot(5,1,2)
plot(timeaxis, Tb_savg_1, '--', 'linewidth', 1.5)
hold on
plot(timeaxis, Tb_savg_2, '--', 'linewidth', 1.5)
plot(timeaxis, Tb_savg_3, '-', 'linewidth', 1.5)
plot(timeaxis, Tb_savg_4, '-.', 'linewidth', 1.5)
plot(timeaxis, Tb_savg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1, str_caption_2, str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim([0,5])
ylabel('T_b(K)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('B. Boundary layer temperature, domain averaged','FontWeight','Bold','FontSize',18)

% Tf
subplot(5,1,3)
plot(timeaxis, Tf_savg_1, '--', 'linewidth', 1.5)
hold on
plot(timeaxis, Tf_savg_2, '--', 'linewidth', 1.5)
plot(timeaxis, Tf_savg_3, '-', 'linewidth', 1.5)
plot(timeaxis, Tf_savg_4, '-.', 'linewidth', 1.5)
plot(timeaxis, Tf_savg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1, str_caption_2, str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim([0,5])
ylabel('T_f(K)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('C. Free troposphere temperature, domain averaged','FontWeight','Bold','FontSize',18)

% shallow cloud
subplot(5,1,4)
plot(timeaxis, bcloud_savg_1, '--', 'linewidth', 1.5)
hold on
plot(timeaxis, bcloud_savg_2, '--', 'linewidth', 1.5)
plot(timeaxis, bcloud_savg_3, '-', 'linewidth', 1.5)
plot(timeaxis, bcloud_savg_4, '-.', 'linewidth', 1.5)
plot(timeaxis, bcloud_savg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1, str_caption_2, str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim([0,5])
ylabel('\sigma_s','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('D. Shallow cloud fraction, domain averaged','FontWeight','Bold','FontSize',18)

% deep cloud 
subplot(5,1,5)
plot(timeaxis, dcloud_savg_1, '--', 'linewidth', 1.5)
hold on
plot(timeaxis, dcloud_savg_2, '--', 'linewidth', 1.5)
plot(timeaxis, dcloud_savg_3, '-', 'linewidth', 1.5)
plot(timeaxis, dcloud_savg_4, '-.', 'linewidth', 1.5)
plot(timeaxis, dcloud_savg_5, '-.', 'linewidth', 1.5)
legend(str_caption_1, str_caption_2, str_caption_3,str_caption_4,str_caption_5,'location','best')
xlim([0,5])
ylabel('\sigma_f','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('E. Deep cloud fraction, domain averaged','FontWeight','Bold','FontSize',18)

savefig(fig_savg_varbs,'[sens][TS]varbs.fig');
close(fig_savg_varbs)

cd ..

%% SOME HELP FUNCTIONS

function [] = write_axis_x()

xticks([0,2000,4000,6000,8000,10000])
xticklabels({'0','2000','4000','6000','8000','10000'})

end
