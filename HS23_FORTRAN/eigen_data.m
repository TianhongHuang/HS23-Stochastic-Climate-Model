%--------------------------------------------------------------------------
% [SBR]: EIGEN_DATA
%--------------------------------------------------------------------------
% author:   Tianhong Huang
% date:     May 2023
% contains: Evec_dat, Eval_ex1, Eval_ex2
% parallel: used in parallel code, # of *.dat = 4 * # of threads
% program:  SH23
% version:  v5.1
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% parameters from model 
%--------------------------------------------------------------------------
dt = input('  [time step] dt = ','s');
dt = str2double(dt);
dx = input('  [spatial size] dx = ','s');
dx = str2double(dx)*1e+3;
xlen = input('  [grid size x] xlen = ','s');
xlen = str2double(xlen);
ylen = input('  [grid size y] ylen = ','s');
ylen = str2double(ylen);
Lx = dx*(xlen/2);
Ly = dx*(ylen/2);
sfx = pi/Lx;
sfy = pi/Ly;
HT = 15500;
hb = 1000;
Cd = 0.025;
Up = 2;
theta_ref = 300;
g = 9.8;
B = 10^(-4)*theta_ref/g;                                              
alpha1 = (HT/pi)*(g/theta_ref);                                            
alpha2 = (HT/pi)*B;                                                        

%--------------------------------------------------------------------------
% define output eigan matrice
%--------------------------------------------------------------------------
k = (0:xlen-1)-xlen/2;
l = (0:ylen-1)-ylen/2;
A = cell(xlen,ylen);
[ V, D, Vinv ]   = deal(cell(xlen,ylen));
[ SIG1, SIG2 ]   = deal(cell(xlen,ylen));
[ diag1, diag2 ] = deal(zeros(1,3));

%--------------------------------------------------------------------------
% eigan value decomposition
%--------------------------------------------------------------------------
for i = 1: xlen
    for j = 1: ylen
        Fmode = [k(i),l(j)];
        if isequal(Fmode,[0,0]) || isequal(Fmode,[0,-ylen/2]) || ...
           isequal(Fmode,[-xlen/2,0]) || isequal(Fmode,[-xlen/2,-ylen/2])
            D{i,j} = zeros(3,3);
            V{i,j} = zeros(3,3);
            Vinv{i,j} = zeros(3,3);
            SIG1{i,j} = zeros(3,3);
            SIG2{i,j} = zeros(3,3);
        else
            A{i,j} = [ 0, -alpha1*((1i*sfx*k(i))^2+(1i*sfy*l(j))^2), 0;...
                 -alpha2, 0, -alpha2*((1i*sfx*k(i))^2+(1i*sfy*l(j))^2)/sqrt(2)/HT;...
                  0, -alpha1*2*sqrt(2)*hb*HT/(hb+HT), Cd*Up*HT/(HT+hb)/hb ];
            [V{i,j},D{i,j}] = eig(A{i,j});
            Vinv{i,j} = V{i,j}\eye(3,3);
            for s = 1: 3
                diag1(s) = (1-exp(-D{i,j}(s,s)*dt))/D{i,j}(s,s);
                diag2(s) = exp(-D{i,j}(s,s)*dt);
            end
            SIG1{i,j} = diag(diag1);
            SIG2{i,j} = diag(diag2);
        end
    end
end


%--------------------------------------------------------------------------
% rescaling eigen pairs & storing eigen data
%--------------------------------------------------------------------------
str_proc_max = input('  [# of proc] what is the max number of processors: ','s');
proc_max = str2double(str_proc_max);
str_proc_max = sprintf('%02.f', proc_max);

cd([pwd,'/eigen_info/'])

for proc_id = 1: proc_max
    
    str_proc_id = num2str(proc_id-1,'%02d');
    
    sub_ylen = ylen/proc_max;
    
    col_ref_idx = (proc_id-1) * sub_ylen;
    
    % define V.dat
    fid = fopen(['V_',str_proc_id,'_',str_proc_max,'.dat'], 'wt');
    for i = 1: xlen
        for j = 1: sub_ylen
            real_part = reshape(real(V{i,j+col_ref_idx}).', ...
                1, numel(V{i,j+col_ref_idx}));
            imag_part = reshape(imag(V{i,j+col_ref_idx}).', ...
                1, numel(V{i,j+col_ref_idx}));
            print_num = [real_part; imag_part];
            fprintf(fid, '(%.17e, %.17e)\n', print_num);
        end
    end
    fclose(fid);
    
    % define Vinv.dat
    fid = fopen(['Vinv_',str_proc_id,'_',str_proc_max,'.dat'], 'wt');
    for i = 1: xlen
        for j = 1: sub_ylen
            real_part = reshape(real(Vinv{i,j+col_ref_idx}).', ...
                1, numel(Vinv{i,j+col_ref_idx}));
            imag_part = reshape(imag(Vinv{i,j+col_ref_idx}).', ...
                1, numel(Vinv{i,j+col_ref_idx}));
            print_num = [real_part; imag_part];
            fprintf(fid, '(%.17e, %.17e)\n', print_num);
        end
    end
    fclose(fid);
    
    % define diag_1.dat
    fid = fopen(['diag1_',str_proc_id,'_',str_proc_max,'.dat'], 'wt');
    for i = 1: xlen
        for j = 1: sub_ylen
            real_part = reshape(real(SIG1{i,j+col_ref_idx}).', ...
                1, numel(SIG1{i,j+col_ref_idx}));
            imag_part = reshape(imag(SIG1{i,j+col_ref_idx}).', ...
                1, numel(SIG1{i,j+col_ref_idx}));
            print_num = [real_part; imag_part];
            fprintf(fid, '(%.17e, %.17e)\n', print_num);
        end
    end
    fclose(fid);
    
    % define diag_2.dat
    fid = fopen(['diag2_',str_proc_id,'_',str_proc_max,'.dat'], 'wt');
    for i = 1: xlen
        for j = 1: sub_ylen
            real_part = reshape(real(SIG2{i,j+col_ref_idx}).', ...
                1, numel(SIG2{i,j+col_ref_idx}));
            imag_part = reshape(imag(SIG2{i,j+col_ref_idx}).', ...
                1, numel(SIG2{i,j+col_ref_idx}));
            print_num = [real_part; imag_part];
            fprintf(fid, '(%.17e, %.17e)\n', print_num);
        end
    end
    fclose(fid);
    
end

cd ..
