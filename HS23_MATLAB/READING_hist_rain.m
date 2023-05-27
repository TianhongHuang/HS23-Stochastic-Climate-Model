%--------------------------------------------------------------------------
% [SBR]: READING_HIST_RAIN
% [TAG]: model variables time-avg/time-meridionally-avg data
% [INTRO]: reading outputs data and storing correspond varbs as .mat file
%          plotting histogram data related w./ rainfall 
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

dtbinnum = param.dtbinnum; % # of bins (dry spell duration)
rtbinnum = param.rtbinnum; % # of bins (rain events duration)
rsbinnum = param.rsbinnum; % # of bins (rain events size)
stbinnum = param.stbinnum; % # of bins (shallow cluster duration)
DryTimeMIN = param.drytimemin; % lower threshold (dry spell duration)
DryTimeMAX = param.drytimemax; % upper threshold (dry spell duration)
RainTimeMIN = param.raintimemin; % lower threshold (rain events duration)
RainTimeMAX = param.raintimemax; % upper threshold (rain events duration)
RainSizeMIN = param.rainsizemin; % lower threshold (rain events size)
RainSizeMAX = param.rainsizemax; % upper threshold (rain events size)
SigcTimeMIN = param.sigctimemin; % lower threshold (shallow cluster duration)
SigcTimeMAX = param.sigctimemax; % upper threshold (shallow cluster duration)
proc_num = param.procnum;
str_proc_num = num2str(proc_num);

%--------------------------------------------------------------------------
% [**]hist: original histogram data
% [**]hist_re: histogram data with evenly scale under log-sense
% [**]binscale: original bin's ending pts vector
% [**]logbinscale: bin's ending pts vector under log-sense
%--------------------------------------------------------------------------
[ rshist, rshist_re, rsbinscale, rslogbinscale ] = read_histogram( ...
    'RainSize_hist.csv', rsbinnum, RainSizeMIN, RainSizeMAX, proc_num );

[ rthist, rthist_re, rtbinscale, rtlogbinscale ] = read_histogram( ...
    'RainTime_hist.csv', rtbinnum, RainTimeMIN, RainTimeMAX, proc_num );

[ dthist, dthist_re, dtbinscale, dtlogbinscale ] = read_histogram( ...
    'DryTime_hist.csv', dtbinnum, DryTimeMIN, DryTimeMAX, proc_num );

[ sthist, sthist_re, stbinscale, stlogbinscale ] = read_histogram( ...
    'SigcTime_hist.csv', stbinnum, SigcTimeMIN, SigcTimeMAX, proc_num );

%--------------------------------------------------------------------------
% FIG_1: histogram of rain events size
%--------------------------------------------------------------------------
fig_rainsize = figure;

line11 = loglog(rsbinscale, rshist_re,'-s','linewidth',1.5);
set(line11, 'MarkerFaceColor', get(line11,'Color'));
title('rainfall size histogram log-log plot','FontWeight','Bold','FontSize',18)
xlim('auto')
xlabel('rainfall event size (mm)','FontWeight','Bold','FontSize',16)
ylabel('rescaled case number','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

cd([pwd,'/output_fig/'])
title_rshist = '[hist]rainsize';
savefig(fig_rainsize,title_rshist);
close(fig_rainsize)

%--------------------------------------------------------------------------
% FIG_2: histogram of rain/dry spell duration
%--------------------------------------------------------------------------
fig_duration = figure;

set(gcf,'units','points','position',[10,10,1200,400])

subplot(1,2,2)
line22 = loglog(dtbinscale, dthist_re,'-s','linewidth',1.5);
set(line22, 'MarkerFaceColor', get(line22,'Color'));
title('dry events duration histogram log-log plot','FontWeight','Bold','FontSize',18)
xlim('auto')
xlabel('dry time (mins)','FontWeight','Bold','FontSize',16)
ylabel('rescaled case number','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

subplot(1,2,1)
line21 = loglog(rtbinscale, rthist_re,'-s','linewidth',1.5);
set(line21, 'MarkerFaceColor', get(line21,'Color'));
title('moist events duration histogram log-log plot','FontWeight','Bold','FontSize',18)
xlim('auto')
xlabel('raining time (mins)','FontWeight','Bold','FontSize',16)
ylabel('rescaled case number','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

title_duration = '[hist]duration';
savefig(fig_duration,title_duration);
close(fig_duration)

%--------------------------------------------------------------------------
% FIG_3: histogram of shallow cluster duration
%--------------------------------------------------------------------------
fig_sigbtime = figure;

line31 = loglog(stbinscale, sthist_re,'-s','linewidth',1.5);
set(line31, 'MarkerFaceColor', get(line31,'Color'));
title('shallow cloud duration histogram log-log plot','FontWeight','Bold','FontSize',18)
xlim('auto')
xlabel('shallow cloud duration (mins)','FontWeight','Bold','FontSize',16)
ylabel('rescaled case number','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

savefig(fig_sigbtime,'[hist]sigbtime');
close(fig_sigbtime)

cd ..

Prs = polyfit(rslogbinscale(1:rsbinnum-20),log(rshist_re(1:rsbinnum-20)),1);
Prt = polyfit(rtlogbinscale(1:rtbinnum-10),log(rthist_re(1:rtbinnum-10)),1);
Pdt = polyfit(dtlogbinscale(1:dtbinnum-10),log(dthist_re(1:dtbinnum-10)),1);
Pst = polyfit(stlogbinscale(1:stbinnum-10),log(sthist_re(1:stbinnum-10)),1);

%--------------------------------------------------------------------------
% save histogram data as 'hist_rain_events.mat'
%--------------------------------------------------------------------------
save('hist_rain_events.mat')

%--------------------------------------------------------------------------
% display important histogram statistics
%--------------------------------------------------------------------------
fprintf('----------------------------------------------------------------------\n')
fprintf('  [histogram][moist-data] figure closed.')
fprintf('\n----------------------------------------------------------------------\n')
fprintf('  [dry time] # of head cases = %.0f.\n', dthist(1))   
fprintf('  [dry time] # of tail cases = %.0f.\n', dthist(dtbinnum))  
fprintf('  [dry time] # of log-log slope = %1.2f.\n', Pdt(1)) 
fprintf('  [rain time] # of head cases = %.0f.\n', rthist(1))   
fprintf('  [rain time] # of tail cases = %.0f.\n', rthist(rtbinnum))  
fprintf('  [rain time] # of log-log slope = %1.2f.\n', Prt(1)) 
fprintf('  [rain size] # of head cases = %.0f.\n', rshist(1))   
fprintf('  [rain size] # of tail cases = %.0f.\n', rshist(rsbinnum))  
fprintf('  [rain size] # of log-log slope = %1.2f.\n', Prs(1))    
fprintf('  [sigc time] # of head cases = %.0f.\n', sthist(1))   
fprintf('  [sigc time] # of tail cases = %.0f.\n', sthist(stbinnum))  
fprintf('  [sigc time] # of log-log slope = %1.2f.\n', Pst(1)) 
fprintf('----------------------------------------------------------------------\n\n')

%--------------------------------------------------------------------------
% SOME HELP FUNCTIONS
%--------------------------------------------------------------------------
function[ hist, hist_rescale, binscale, logbinscale ] = read_histogram( ...
    filename, binnum, binmin, binmax, proc_num )

cd([pwd,'/output_parallel'])

str_proc_num = num2str(proc_num);

for proc_id = 1: proc_num
    str_proc_id = num2str(proc_id-1,'%02d');
    str_f90_op  = ['/MAX_',str_proc_num,'_ID_',str_proc_id,'/stats/'];
    cd([pwd,str_f90_op])
    if proc_id == 1
        hist_init = csvread(filename);
        hist = hist_init(:,1:size(hist_init,2)-1);
    else
        hist_temp = csvread(filename);
        hist_temp = hist_temp(:,1:size(hist_temp,2)-1);
        hist = hist + hist_temp;
    end
    cd ..
    cd ..
end

cd ..

logbinsize = (log(binmax)-log(binmin))/(binnum-1);
logbinscale = (0:binnum-2)*logbinsize+log(binmin);
binscale = exp(logbinscale);
binscale_shift = circshift(binscale, [0,-1]);
binscale_shift(end) = binmax;
binwidth = binscale_shift - binscale;
hist_rescale = hist(1:binnum-1)./binwidth';

end
