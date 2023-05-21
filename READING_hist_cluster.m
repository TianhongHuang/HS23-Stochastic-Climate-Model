%--------------------------------------------------------------------------
% [SBR]: READING_HIST_CLUSTER
% [TAG]: model variables time-avg/time-meridionally-avg data
% [INTRO]: reading outputs data and storing correspond varbs as .mat file
%          plotting histogram data related w./ cloud clusters
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

proc_num = param.procnum;
xlen = param.xlen;
ylen = param.ylen;

str_proc_max = num2str(proc_num);

data_time = csvread([pwd,'/output_parallel/time_sequence.csv']);
data_time = data_time(:,1)';
data_size = numel(data_time);

%--------------------------------------------------------------------------
% init_setup for histogram data
%--------------------------------------------------------------------------
max_size = xlen*ylen;
min_size = 1;
binnum = 30;
logbinsize = (log(max_size)-log(min_size))/(binnum-1);
logbinscale = log(min_size):logbinsize:log(max_size);
binscale = exp(logbinscale);
binscale_shift = circshift(binscale(1:binnum-1), [0,-1]);
binscale_shift(end) = max_size;
binwidth = binscale_shift - binscale(1:binnum-1);
csbinscale = binscale(1:binnum-1);

scs_hist = zeros(1, binnum); % shallow-cluster-size histogram data
dcs_hist = zeros(1, binnum); % deep-convective-size histogram data

%--------------------------------------------------------------------------
% collecting cluster size histogram from FORTRAN outputs
%--------------------------------------------------------------------------
for t = 1: data_size
    
    str_step_num = num2str(t-1, '%04.f');
    
    shallow_cloud = read_parallel_data('shallow-cloud-ast', str_proc_max, str_step_num);
    deep_cloud = read_parallel_data('deep-cloud', str_proc_max, str_step_num);
    shallow_clusters = Cluster_Count_v2(shallow_cloud);
    deep_clusters = Cluster_Count_v2(deep_cloud);
    num_of_clouds_1 = numel(shallow_clusters);
    num_of_clouds_2 = numel(deep_clusters);
    
    for i = 1: num_of_clouds_1
        index1 = floor((log(shallow_clusters(i))-log(min_size))/logbinsize)+1;
        index1 = min(index1, binnum);
        scs_hist(index1) = scs_hist(index1)+1;
    end
    
    for i = 1: num_of_clouds_2
        index2 = floor((log(deep_clusters(i))-log(min_size))/logbinsize)+1;
        index2 = min(index2, binnum);
        dcs_hist(index2) = dcs_hist(index2)+1;
    end
    
    if mod(t-1,10) == 0
        fprintf('  cluster data scanned: %2.2f %%.\n', ((t-1)/data_size)*1E+2)
    end
    
end

scs_hist_rescale = scs_hist(1:binnum-1)./binwidth;
dcs_hist_rescale = dcs_hist(1:binnum-1)./binwidth;

%--------------------------------------------------------------------------
% FIG: cluster size histogram data
%--------------------------------------------------------------------------
fig_clustersize = figure;

set(gcf,'units','points','position',[10,10,1200,400])

subplot(1,2,1)
line1 = loglog(csbinscale(2:end),scs_hist_rescale(2:end),'-s','linewidth',1.5);
set(line1, 'MarkerFaceColor', get(line1,'Color'));
title('shallow cluster size histogram log-log plot','FontWeight','Bold','FontSize',18)
xlim('auto')
xlabel('cluster size (num of gird pts)','FontWeight','Bold','FontSize',16)
ylabel('rescaled case number','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

subplot(1,2,2)
line2 = loglog(csbinscale(2:end),dcs_hist_rescale(2:end),'-s','linewidth',1.5);
set(line2, 'MarkerFaceColor', get(line2,'Color'));
title('deep cluster size histogram log-log plot','FontWeight','Bold','FontSize',18)
xlim('auto')
xlabel('cluster size (num of gird pts)','FontWeight','Bold','FontSize',16)
ylabel('rescaled case number','FontWeight','Bold','FontSize',16)
ax = gca;
ax.FontSize = 16;

cd([pwd,'/output_fig/'])
title_cshist = '[hist]cluster';
savefig(fig_clustersize,title_cshist);
close(fig_clustersize)
cd ..

Pscs = polyfit(logbinscale(2:binnum-5),log(scs_hist_rescale(2:binnum-5)),1);
Pdcs = polyfit(logbinscale(2:binnum-5),log(dcs_hist_rescale(2:binnum-5)),1);

fprintf('  [histogram][cluster-size] figure closed.\n')
fprintf('  [shallow][cluster size] log-log slope = %1.2f.\n', Pscs(1)) 
fprintf('  [deep][cluster size] log-log slope = %1.2f.\n', Pdcs(1)) 

%--------------------------------------------------------------------------
% save histogram data as 'hist_cluster_size.mat'
%--------------------------------------------------------------------------
save('hist_cluster_size.mat')

fprintf('----------------------------------------------------------------------\n\n')

%--------------------------------------------------------------------------
% HELP FUNCTION: reading parallel data from FORTRAN programming
%--------------------------------------------------------------------------
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

%--------------------------------------------------------------------------
% HELP FUNCTION: counting cluster size from outputs data
%--------------------------------------------------------------------------
function [Cloud_Counts] = Cluster_Count_v2(q)

% Compute total number of clouds
% This is the largest possible cloud in the realization

[N,M] = size(q);
% N = num of rows
% M = num of cols
I = find(q>0);
% I = indice of all clouds
% Num_clouds = length(I);

tic

buffer = [];
% buffer = pocket of curr cloud grid pts
neighbor = zeros(4,1);
% neighbor = N/W/E/S grid pts of curr cloud 
% IMPORTANT: in modifed version, we DO NOT count neighbor of boundary pts
Clouds_in_Cluster = [];
% C_i_C = array that counts every single cloud clusters
% C_i_C(i) = num of grid pts (size) of the i-th cloud cluster
i = 0;

while isempty(I)==0
   
    i = i+1;

    curr_cloud = I(1);

    buffer = curr_cloud;
    
    % init_size of new cluster = zero
    Clouds_in_Cluster(i) = 0;
    
    while isempty(buffer) == 0
       
        % pick last pt in curr pocket as new center pt
        curr_cloud = buffer(end);

        % Assign the neighbor corrdinates
        neighbor(1) = curr_cloud+1; % South
        neighbor(2) = curr_cloud-1; % North
        neighbor(3) = curr_cloud-N; % West
        neighbor(4) = curr_cloud+N; % East
    
        for j=1:4
            % DO NOT count neighbor pts that out of domain
            if neighbor(j) < 1 || neighbor(j) > M*N
                flag_bdry = 1;
            % DO NOT count south neighbor of bottom side pts  
            elseif mod(curr_cloud,N) == 0 && j == 1
                flag_bdry = 1;
            % DO NOT count north neighbor of top side pts
            elseif mod(curr_cloud-1,N) == 0 && j == 2
                flag_bdry = 1;
            % green flag for all other cases    
            else
                flag_bdry = 0;
            end
            
            % Checking conditions for green flag case
            if flag_bdry == 0                
                
                if q(neighbor(j))>0 ...
                        && isempty(buffer(buffer==neighbor(j))) ...
                        && isempty(I(I==neighbor(j)))==0 ...
                    % cond_1: q(neighbor(j))>0 indicates neighbor cloud
                    % cond_2: isempty(buffer(buffer==neighbor(j)))
                    % checking cloud pt was not added into pocket before
                    % cond_3: isempty(I(I==neighbor(j)))==0
                    % checking cloud pt was not counted into former clusters
                    
                    % if q(neighbor(j)) s.t cond_1-3 then add into pocket
                    buffer = [buffer, neighbor(j)];
                    
                end
            end
        end
       
        % Increment counter for cluster
        Clouds_in_Cluster(i) = Clouds_in_Cluster(i)+1;
        
        % Boot cloud out of buffer
        buffer = buffer(buffer~=curr_cloud);
        
        % Boot cloud from Indices
        I = I(I~=curr_cloud);
        
    end
    
    
end

Cloud_Counts = Clouds_in_Cluster;

end
