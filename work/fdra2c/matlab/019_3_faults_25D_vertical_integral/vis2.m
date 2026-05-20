% Yudong May 7 2025
% load the data for the case twosigmas
addpath(genpath('/home/yudong/FDRA2C_external/fdra2c/'))
%addpath(genpath('/home/yudong/FDRA2C/matlab-misc/'))
addpath(genpath('/home/yudong/FDRA2C_external_GF/fdra2c/work/fdra2c/matlab/'))
data_dir = '/nfs/camcat001/yudong/fdra2c/external_GF/test19/';
sim_name = 'v2Gn';
load([data_dir sim_name '.mat'])     
q = fdra2c('qload', s.fs.cf,1);
% Yudong May 22 2025, load the stress 
q2 = fdra2c('qload2', s,1);
%snew = replace_path_struct(s, s.fs.cf.dir, data_dir)
%q = fdra2c('qload', snew.fs.cf,1);  
% q = q2load([data_dir sim_name],1);
% fdra2c('vi_Start', q);  

