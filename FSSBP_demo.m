% A demo for BP_I, BP_T, SSBP, FBP and FSSBP.
% FSSBP: Fast Spatial-Spectral Back Projection Based on Pan-sharpening
% Iterative Optimization, 2023. https://doi.org/10.3390/rs15184543.

clc;
clear;
close all;
file_path = matlab.desktop.editor.getActive;
cd(fileparts(file_path.Filename));

n_run_times = 1;
update_output_qi = 1;
data_name = 'D-none_16';%'F-QB_1' 'F-GeoEye1_12' 'D-none_16' 'D-WV2_19'

% In order to simplify the demonstration process, all mat files are
% accompanied by an initial solution to the base method as well as a number
% of associated global parameters, in addition to the data itself.

load(['Datasets/' data_name '.mat']);
switch(data_name)
    case 'D-none_16'
        base_method_name = 'ATWT-M3';
        I_F0 = I_F_ATWT_M3;
    case 'D-WV2_19'
        base_method_name = 'PNN';
        I_F0 = I_F_PNN;
    case 'F-QB_1'
        base_method_name = 'MTF-GLP';
        I_F0 = I_F_MTF_GLP;
    case 'F-GeoEye1_12'
        base_method_name = 'GFPCA';
        I_F0 = I_F_GFPCA;
end

n_methods = 0;
BP_case_name = {'BP_T', 'BP_I', 'SSBP', 'FBP', 'FSSBP'};
for ii = 1:length(BP_case_name)
    BP_Opts = init_BP_options(BP_case_name{ii}, is_FS);
    
    n_methods = n_methods+1;
    Method_{n_methods} = [BP_Opts.name '_' base_method_name ...
        '_RT_' num2str(BP_Opts.RT_mode) '_gamma' ...
        num2str(BP_Opts.gamma), '_tau', ...
        num2str(BP_Opts.m_lambda) '_DT', ...
        num2str(BP_Opts.DT_mode)];
    t2=tic;    
    for i = 1:n_run_times
        I_F_BP{n_methods} = BP_Wrapper(double(I_F0), I_MS_LR, I_PAN, ratio, sensorInf, BP_Opts);
    end
    Time_{n_methods}=toc(t2)/n_run_times;
    fprintf('Elaboration time %s: %.4f [sec]\n',Method_{n_methods},Time_{n_methods});
end

dat_time = datestr(now,30);

if is_FS
    I_GT = [];
end

QI_ = cell(1, n_methods);
for i = 1:n_methods
    QI_{i} = indices_eval_wrapper(I_F_BP{i}, I_GT, I_MS, I_MS_LR, I_PAN, ratio,...
        L, sensorInf, is_FS, Qblocks_size, flag_cut_bounds, dim_cut,thvalues);    
end
QI_matrix = QI2matrix(QI_,Method_,Time_,is_FS);
if update_output_qi
    disp( ['Writting xls of Dataset: ' data_name]);
    xlswrite(['Results/qi/' datestr(now,30) '_' data_name], QI_matrix, ['Dataset ' data_name]);
end
