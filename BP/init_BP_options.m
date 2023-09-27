function BP_Opts = init_BP_options(type, is_FS)
BP_Opts.name = type;

if is_FS
    BP_Opts.DT_mode = 1;
    BP_Opts.gamma = 16;
    BP_Opts.RT_mode = 1;
    BP_Opts.m_lambda = 1;
else
    BP_Opts.DT_mode = 0;
    BP_Opts.gamma = 1;
    BP_Opts.RT_mode = 2;
    BP_Opts.m_lambda = 0.1;
end

BP_Opts.n_reg = 100;
%BP_Opts.m_lambda = 1;
BP_Opts.auto_m_lambda = 0;
switch (type)
    case 'BP_I'
        BP_Opts.is_iter = 1;
        BP_Opts.is_SS = 0;
        BP_Opts.gamma = 1; % Step size fixed at 1.
        BP_Opts.DT_mode = 0; % Use of ideal interpolation.
    case 'BP_T'
        BP_Opts.is_iter = 1;
        BP_Opts.is_SS = 0;
        BP_Opts.gamma = 1; % Step size fixed at 1.
        BP_Opts.DT_mode = 1; % Use of degraded transposition.
    case 'FBP'
        BP_Opts.is_iter = 0;
        BP_Opts.is_SS = 0;
    case 'SSBP'
        BP_Opts.is_iter = 1;
        BP_Opts.is_SS = 1;
    case 'ASSBP'
        BP_Opts.is_iter = 1;
        BP_Opts.is_SS = 1;
        BP_Opts.auto_m_lambda = 1; 
    case 'FSSBP'
        BP_Opts.is_iter = 0;
        BP_Opts.is_SS = 1;
    case 'AFSSBP'
        BP_Opts.is_iter = 0;
        BP_Opts.is_SS = 1;
        BP_Opts.auto_m_lambda = 1;
end

if ~is_FS
    BP_Opts.mu = 0.0098;
else
    BP_Opts.mu = 0.2;
end

