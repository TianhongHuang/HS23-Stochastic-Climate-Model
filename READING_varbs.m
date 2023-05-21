%--------------------------------------------------------------------------
% [SBR]: READING_VARBS
% [TAG]: model variables time-avg/time-meridionally-avg data
% [INTRO]: reading outputs data and storing correspond varbs as .mat file
%          plotting time-averaged model variables
%          plotting space-and-time evolution of model variables
%--------------------------------------------------------------------------
% @author:  Tianhong Huang
% @date:    May 2023
% @program: SH23
% @version: v5.1
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% using 'read_namelist.m' to input NAMELIST file as parameters
%--------------------------------------------------------------------------
%   Written by:     S.Lazerson (lazerson@pppl.gov)
%   Version:        1.5
%   Date:           9/13/12
%--------------------------------------------------------------------------
param = read_namelist('NAMELIST','model');

dt = param.dt;
proc_num = param.procnum;
domain_x = param.ylen;
domain_y = param.xlen;
str_proc_max = num2str(proc_num);
data_time = csvread([pwd,'/output_parallel/time_sequence.csv']);
data_time = data_time(:,1)';
data_size = numel(data_time);

mkdir output_fig % save plot figures in folder [output_fig]

%% READING TIME-AVGED DATA & PLOTTING 

% ocean temperature
To_tavg = read_parallel_data('stats', str_proc_max, 'To_tavg');
% boundary layer temperature
Tb_tavg = read_parallel_data('stats', str_proc_max, 'Tb_tavg');
% free troposphere temperature
Tf_tavg = read_parallel_data('stats', str_proc_max, 'Tf_tavg');
% baroclinic mode velocity
u1_tavg = read_parallel_data('stats', str_proc_max, 'u1_tavg');
v1_tavg = read_parallel_data('stats', str_proc_max, 'v1_tavg');
% column water vapors
qf_tavg = read_parallel_data('stats', str_proc_max, 'qf_tavg');
qb_tavg = read_parallel_data('stats', str_proc_max, 'qtb_tavg');
cwv_tavg = qf_tavg + qb_tavg;
% cloud fractions (time-avged)
scloud_rate = read_parallel_data('stats', str_proc_max, 'shallow_cloud_rate');
scloud_tavg = read_parallel_data('stats', str_proc_max, 'shallow_cloud_tavg');
dcloud_rate = read_parallel_data('stats', str_proc_max, 'deep_cloud_rate');
dcloud_tavg = read_parallel_data('stats', str_proc_max, 'deep_cloud_tavg');
bcloud_rate = read_parallel_data('stats', str_proc_max, 'shallow_cloud_ast_rate');
bcloud_tavg = read_parallel_data('stats', str_proc_max, 'shallow_cloud_ast_tavg');
% cloud fractions (time-series)
scloud_rate = sum(scloud_rate,2)/proc_num;
dcloud_rate = sum(dcloud_rate,2)/proc_num;
bcloud_rate = sum(bcloud_rate,2)/proc_num;

%--------------------------------------------------------------------------
% FIG_1.1: time averaged model variables plot
%--------------------------------------------------------------------------
fig_tavg = figure;
set(gcf,'units','points','position',[10,10,1000,1600])

subplot(8,1,1)
pcolor(To_tavg);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('A. T_o(K)','FontWeight','Bold','FontSize',18)

subplot(8,1,2)
pcolor(Tb_tavg);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('B. T_b(K)','FontWeight','Bold','FontSize',13)

subplot(8,1,3)
pcolor(qb_tavg*1e+3);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('C. q_b(mm)','FontWeight','Bold','FontSize',13)

subplot(8,1,4)
pcolor(bcloud_tavg);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('D. \sigma_b','FontWeight','Bold','FontSize',13)

subplot(8,1,5)
pcolor(Tf_tavg);
shading interp
Tfclim = mean(mean(Tf_tavg));
caxis([Tfclim-0.1 Tfclim+0.1])
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('E. T_f(K)','FontWeight','Bold','FontSize',13)

subplot(8,1,6)
pcolor(qf_tavg*1e+3);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('F. q_f(mm)','FontWeight','Bold','FontSize',13)

subplot(8,1,7)
pcolor(dcloud_tavg);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
title('G. \sigma_f','FontWeight','Bold','FontSize',13)

subplot(8,1,8)
pcolor(v1_tavg);
shading interp
colorbar
write_axis_x()
write_axis_y()
ax = gca;
ax.FontSize = 16;
xlabel('x (km)','FontSize',16)
title('H. u_1(m/s)','FontWeight','Bold','FontSize',13)

cd([pwd,'/output_fig/'])
savefig(fig_tavg,'[TA]varbs');
close(fig_tavg)
cd ..
fprintf('----------------------------------------------------------------------\n')
fprintf('  [TA][varbs] figure closed.\n')

%% TIME EVOLUTION DATA & PLOTTING
[xlen,ylen] = size(To_tavg);
To_yavg = zeros(data_size,ylen);
Tb_yavg = zeros(data_size,ylen);
Tf_yavg = zeros(data_size,ylen);
qb_yavg = zeros(data_size,ylen);
qf_yavg = zeros(data_size,ylen);
v1_yavg = zeros(data_size,ylen);
cwv_yavg = zeros(data_size,ylen);
scloud_yavg = zeros(data_size,ylen);
dcloud_yavg = zeros(data_size,ylen);
bcloud_yavg = zeros(data_size,ylen);

for t = 1: data_size
    
    str_step_num = num2str(t-1, '%04.f');
    To = read_parallel_data('To', str_proc_max, str_step_num);
    Tb = read_parallel_data('Tb', str_proc_max, str_step_num);
    Tf = read_parallel_data('Tf', str_proc_max, str_step_num);
    v1 = read_parallel_data('v1', str_proc_max, str_step_num);
    qf = read_parallel_data('qf', str_proc_max, str_step_num);
    qtb = read_parallel_data('qtb', str_proc_max, str_step_num);
    scloud = read_parallel_data('shallow-cloud', str_proc_max, str_step_num);
    dcloud = read_parallel_data('deep-cloud', str_proc_max, str_step_num);
    bcloud = read_parallel_data('shallow-cloud-ast', str_proc_max, str_step_num);
    
    for j = 1: ylen
        To_yavg(t,j) = mean(To(:,j));
        Tb_yavg(t,j) = mean(Tb(:,j));
        Tf_yavg(t,j) = mean(Tf(:,j));
        qb_yavg(t,j) = mean(qtb(:,j));
        qf_yavg(t,j) = mean(qf(:,j));
        v1_yavg(t,j) = mean(v1(:,j));
        cwv_yavg(t,j) = mean(qtb(:,j)+qf(:,j));
        scloud_yavg(t,j) = mean(scloud(:,j));
        dcloud_yavg(t,j) = mean(dcloud(:,j));
        bcloud_yavg(t,j) = mean(bcloud(:,j));
    end
    
end

%--------------------------------------------------------------------------
% saving time-avg & time-meridionally-avg data as .mat file
%--------------------------------------------------------------------------
save('stats_varbs_ast.mat','To_tavg','Tb_tavg','Tf_tavg','u1_tavg','v1_tavg',...
    'qf_tavg','qb_tavg','scloud_tavg','dcloud_tavg','scloud_rate',...
    'dcloud_rate','bcloud_tavg','bcloud_rate','To_yavg','Tb_yavg',...
    'Tf_yavg','qf_yavg','qb_yavg','v1_yavg','cwv_yavg','scloud_yavg',...
    'dcloud_yavg','bcloud_yavg','data_time')

fprintf('  [mat][varbs] outputs done.\n')

%--------------------------------------------------------------------------
% FIG_2.1: time evolution of model variables
%--------------------------------------------------------------------------
data_size_sub = data_size;
data_size = 0.99*data_size;
stats_tic = 0.5;
time_label = [0,round(0.2*round((1-stats_tic)*data_size)),...
    2*round(0.2*round((1-stats_tic)*data_size)),...
    3*round(0.2*round((1-stats_tic)*data_size)),...
    4*round(0.2*round((1-stats_tic)*data_size))];
time_index = [round(stats_tic*data_size)+1,...
    round(stats_tic*data_size)+round(0.2*round((1-stats_tic)*data_size)),...
    round(stats_tic*data_size)+2*round(0.2*round((1-stats_tic)*data_size)),...
    round(stats_tic*data_size)+3*round(0.2*round((1-stats_tic)*data_size)),...
    round(stats_tic*data_size)+4*round(0.2*round((1-stats_tic)*data_size))];

fig_varbs_evol = figure;
set(gcf,'units','points','position',[10,10,1000,1600])

subplot(6,1,1)
pcolor(To_yavg(round(stats_tic*data_size)+1:data_size,:));
shading interp
colorbar
write_axis_x()
yticks(time_label)
yticklabels({'0',num2str(data_time(time_index(2))/365,'%4.1f'),...
    num2str(data_time(time_index(3))/365,'%4.1f'),...
    num2str(data_time(time_index(4))/365,'%4.1f'),...
    num2str(data_time(time_index(5))/365,'%4.1f')})
ylabel('time (years)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('A. meridional avg T_o(K)','FontWeight','Bold','FontSize',18)

subplot(6,1,2)
pcolor(qb_yavg(round(stats_tic*data_size)+1:data_size,:)*1e+3);
shading interp
colorbar
write_axis_x()
yticks(time_label)
yticklabels({'0',num2str(data_time(time_index(2))/365,'%4.1f'),...
    num2str(data_time(time_index(3))/365,'%4.1f'),...
    num2str(data_time(time_index(4))/365,'%4.1f'),...
    num2str(data_time(time_index(5))/365,'%4.1f')})
ylabel('time (years)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('B. meridional avg q_b(mm)','FontWeight','Bold','FontSize',18)

subplot(6,1,3)
pcolor(bcloud_yavg(round(stats_tic*data_size)+1:data_size,:));
shading interp
colorbar
write_axis_x()
yticks(time_label)
yticklabels({'0',num2str(data_time(time_index(2))/365,'%4.1f'),...
    num2str(data_time(time_index(3))/365,'%4.1f'),...
    num2str(data_time(time_index(4))/365,'%4.1f'),...
    num2str(data_time(time_index(5))/365,'%4.1f')})
ylabel('time (years)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('C. meridional avg \sigma_b','FontWeight','Bold','FontSize',18)

subplot(6,1,4)
pcolor(qf_yavg(round(stats_tic*data_size)+1:data_size,:)*1e+3);
shading interp
colorbar
write_axis_x()
yticks(time_label)
yticklabels({'0',num2str(data_time(time_index(2))/365,'%4.1f'),...
    num2str(data_time(time_index(3))/365,'%4.1f'),...
    num2str(data_time(time_index(4))/365,'%4.1f'),...
    num2str(data_time(time_index(5))/365,'%4.1f')})
ylabel('time (years)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('D. meridional avg q_f(mm)','FontWeight','Bold','FontSize',18)

subplot(6,1,5)
pcolor(dcloud_yavg(round(stats_tic*data_size)+1:data_size,:));
shading interp
colorbar
write_axis_x()
yticks(time_label)
yticklabels({'0',num2str(data_time(time_index(2))/365,'%4.1f'),...
    num2str(data_time(time_index(3))/365,'%4.1f'),...
    num2str(data_time(time_index(4))/365,'%4.1f'),...
    num2str(data_time(time_index(5))/365,'%4.1f')})
ylabel('time (years)','FontSize',16)
ax = gca;
ax.FontSize = 16;
title('E. meridional avg \sigma_f','FontWeight','Bold','FontSize',18)

subplot(6,1,6)
pcolor(v1_yavg(round(stats_tic*data_size)+1:data_size,:));
shading interp
colorbar
write_axis_x()
yticks(time_label)
yticklabels({'0',num2str(data_time(time_index(2))/365,'%4.1f'),...
    num2str(data_time(time_index(3))/365,'%4.1f'),...
    num2str(data_time(time_index(4))/365,'%4.1f'),...
    num2str(data_time(time_index(5))/365,'%4.1f')})
ylabel('time (years)','FontSize',16)
ax = gca;
ax.FontSize = 16;
xlabel('x (km)','FontSize',16)
title('F. meridional avg u_1(m/s)','FontWeight','Bold','FontSize',18)

cd([pwd,'/output_fig/'])
savefig(fig_varbs_evol,'[TE]varbs');
close(fig_varbs_evol)
cd ..

data_size = data_size_sub;
fprintf('  [TE][varbs] figure closed.\n')
fprintf('----------------------------------------------------------------------\n\n')

%% SOME HELP FUNCTIONS

function [output] = read_parallel_data(str_variable, str_proc_max, str_step_num)

proc_max = str2double(str_proc_max);

str_proc_max = num2str(proc_max, '%02d');

if proc_max == 1
    str_f90_op = ['/output_serial/',str_variable];
    cd([pwd,str_f90_op])
    csv_files = dir(pwd);
    for idx = 1: size(csv_files,1)
        if contains(csv_files(idx).name, str_step_num)
            data = csvread(csv_files(idx).name);
            yLen = size(data,2)-1;
            data = data(:,1:yLen);
            workspace_name = [str_variable,'_serial'];
            assignin('base', workspace_name, data)
        end
    end   
    cd ..
    cd ..
    output = data;
end

if proc_max > 1
    cd([pwd,'/output_parallel'])
    for proc_id = 1: proc_max
        str_proc_id = num2str(proc_id-1,'%02d');
        str_f90_op  = ['/MAX_',str_proc_max,'_ID_',...
            str_proc_id,'/',str_variable,'/'];
        cd([pwd,str_f90_op])
        csv_files = dir(pwd);
        if proc_id == 1
            data_global = read_init_thread(csv_files, str_step_num);
        else
            data_global = combine_subgrids(data_global, csv_files, str_step_num);
        end
        cd ..
        cd ..
        clear csv_files
    end
    output = data_global;
    cd ..
end

end

function [output] = read_init_thread( file_dir, step_num )
for idx = 1: size(file_dir,1)
    if contains(file_dir(idx).name, step_num)
        data = csvread(file_dir(idx).name);
        yLen = size(data,2)-1;
        data = data(:,1:yLen);
        output = data;
    end
end
end

function [output] = combine_subgrids( data_global, file_dir, step_num )
for idx = 1: size(file_dir,1)
    if contains(file_dir(idx).name, step_num)
        data = csvread(file_dir(idx).name);
        yLen = size(data,2)-1;
        data = data(:,1:yLen);
        output = [data_global,data];
    end
end
end


function [] = write_axis_x()
    xticks([200,400,600,800,1000,1200,1400,1600,1800,2000])
    xticklabels({'1000','2000','3000','4000','5000','6000','7000','8000','9000','10000'})
    %xlabel('x (km)','FontSize',16)
end

function [] = write_axis_y()
    yticks([0,50,100,150,200])
    yticklabels({'-500','-250','0','250','500'})
    ylabel('y (km)','FontSize',16)
end
