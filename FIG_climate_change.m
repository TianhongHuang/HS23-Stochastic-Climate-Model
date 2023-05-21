%--------------------------------------------------------------------------
% [SBR]: CLIMATE_CHANGE
% [TAG]: climate-change studies
% [INTRO]: comparing of time-meridionally-averaged model variables between
%          standard simulation v.s climate-change simulation
%--------------------------------------------------------------------------
% @author:  Tianhong Huang
% @date:    May 2023
% @program: SH23
% @version: v5.1
%--------------------------------------------------------------------------


%% LOADING OUTPUTS FROM STANDARD SIMULATION
%--------------------------------------------------------------------------
% please name your standard simulations outputs w./ 'stats_varbs_11.mat'
% which can be induced by 'READING-DATA.m'
%--------------------------------------------------------------------------
load('stats_varbs_11.mat')
str_caption_1 = 'standard simulation';

% TIME AVERAGED DATA
To_tavg_1 = To_tavg; % ocean temperature
Tb_tavg_1 = Tb_tavg; % boundary layer temperature
Tf_tavg_1 = Tf_tavg; % free-tropospheric temperature
qb_tavg_1 = qb_tavg; % boundary layer column water vapor
qf_tavg_1 = qf_tavg; % free-tropospheric column water vapor
v1_tavg_1 = v1_tavg; % baroclinic mode of zonal wind velocity
scloud_tavg_1 = scloud_tavg; % shllow cloud fraction
dcloud_tavg_1 = dcloud_tavg; % deep convective cloud fraction
bcloud_tavg_1 = bcloud_tavg; % modified shallow cloud fraction

% TIME-MERIDIONALLY AVERAGED DATA
v1_yavg_1 = mean(v1_tavg_1);
To_yavg_1 = mean(To_tavg_1);
Tb_yavg_1 = mean(Tb_tavg_1);
Tf_yavg_1 = mean(Tf_tavg_1);
qf_yavg_1 = mean(qf_tavg_1);
qb_yavg_1 = mean(qb_tavg_1);
sigc_yavg_1 = mean(scloud_tavg_1);
sigf_yavg_1 = mean(dcloud_tavg_1);
sigb_yavg_1 = mean(bcloud_tavg_1);

%% LOADING OPS FROM CLIMATE-CHANGE SIMULATION
%--------------------------------------------------------------------------
% please name your climate-change simulations outputs w./ 
% 'stats_varbs_12.mat', which can be induced by 'READING-DATA.m'
%--------------------------------------------------------------------------
load('stats_varbs_12.mat')
str_caption_2 = 'climate-change simulation';

% TIME AVERAGED DATA
To_tavg_2 = To_tavg;
Tb_tavg_2 = Tb_tavg;
Tf_tavg_2 = Tf_tavg;
v1_tavg_2 = v1_tavg;
qf_tavg_2 = qf_tavg;
qb_tavg_2 = qb_tavg;
scloud_tavg_2 = scloud_tavg;
dcloud_tavg_2 = dcloud_tavg;
bcloud_tavg_2 = bcloud_tavg;

% TIME-MERIDIONALLY AVERAGED DATA
v1_yavg_2 = mean(v1_tavg_2);
To_yavg_2 = mean(To_tavg_2);
Tb_yavg_2 = mean(Tb_tavg_2);
Tf_yavg_2 = mean(Tf_tavg_2);
qf_yavg_2 = mean(qf_tavg_2);
qb_yavg_2 = mean(qb_tavg_2);
sigc_yavg_2 = mean(scloud_tavg_2);
sigf_yavg_2 = mean(dcloud_tavg_2);
sigb_yavg_2 = mean(bcloud_tavg_2);

%% PLOTTING
%--------------------------------------------------------------------------
% save plot figures in folder [output_fig]
%--------------------------------------------------------------------------
cd([pwd,'/output_fig/'])
xaxis = 5: 5: 10000;

%--------------------------------------------------------------------------
% FIG_1: comparison of time-moridionally-avged model variables
%--------------------------------------------------------------------------
fig_yavg_varbs = figure;
set(gcf,'units','points','position',[10,10,1200,1400])

% meridional avg To
subplot(8,1,1)
plot(xaxis, To_yavg_1, 'linewidth', 1.5)
hold on
plot(xaxis, To_yavg_2, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('T_o(K)','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('A. Ocean temperature','Fontweight','Bold','Fontsize',18)

% meridional avg Tb
subplot(8,1,2)
plot(xaxis, Tb_yavg_1, 'linewidth', 1.5)
hold on
plot(xaxis, Tb_yavg_2, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('T_b(K)','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('B. Boundary layer temperature','Fontweight','Bold','Fontsize',18)

% meridional avg qb
subplot(8,1,3)
plot(xaxis, qb_yavg_1*1e+3, 'linewidth', 1.5)
hold on
plot(xaxis, qb_yavg_2*1e+3, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('q_b(mm)','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('C. Boundary layer column water vapor','Fontweight','Bold','Fontsize',18)

% meridional avg shallow cloud
subplot(8,1,4)
plot(xaxis, sigb_yavg_1, 'linewidth', 1.5)
hold on
plot(xaxis, sigb_yavg_2, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('\sigma_s','FontSize',16)
ylim([0.25 1])
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('D. modified shallow cloud fraction','Fontweight','Bold','Fontsize',18)

% meridional avg Tf
subplot(8,1,5)
plot(xaxis, Tf_yavg_1, 'linewidth', 1.5)
hold on
plot(xaxis, Tf_yavg_2, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('T_f(K)','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('E. Free-tropospheric temperature','Fontweight','Bold','Fontsize',18)

% meridional avg qf
subplot(8,1,6)
plot(xaxis, (qf_yavg_1)*1e+3, 'linewidth', 1.5)
hold on
plot(xaxis, (qf_yavg_2)*1e+3, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('q_f(mm)','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('F. Free-tropospheric column water vapor','Fontweight','Bold','Fontsize',18)

% meridional avg deep cloud
subplot(8,1,7)
plot(xaxis, sigf_yavg_1, 'linewidth', 1.5)
hold on
plot(xaxis, sigf_yavg_2, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('\sigma_f','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
title('G. Deep convective cloud fraction','Fontweight','Bold','Fontsize',18)

% meridional avg zonal wind velocity
subplot(8,1,8)
plot(xaxis, v1_yavg_1, 'linewidth', 1.5)
hold on
plot(xaxis, v1_yavg_2, 'linewidth', 1.5)
legend(str_caption_1,str_caption_2,'location','best')
ylabel('u_1 (m/s)','FontSize',16)
ax = gca;
ax.FontSize = 16;
write_axis_x()
xlabel('x (km)','FontSize',16)
title('H. First baroclinic mode of zonal wind velocity','Fontweight','Bold','Fontsize',18)

savefig(fig_yavg_varbs,'[climate-change][TMA]varbs.fig');
close(fig_yavg_varbs)

cd ..

%% SOME HELP FUNCTIONS

function [] = write_axis_x()

xticks([0,2000,4000,6000,8000,10000])
xticklabels({'0','2000','4000','6000','8000','10000'})

end
