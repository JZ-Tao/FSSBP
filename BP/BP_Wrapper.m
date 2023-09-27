% FSSBP: Fast Spatial-Spectral Back Projection Based on Pan-sharpening
% Iterative Optimization, 2023. https://doi.org/10.3390/rs15184543.
% Inputs: 
%   I_0: Initial sharpening result.
%   I_MS_LR: LR MS image.
%   I_PAN: HR PAN image.
%   ratio: The resolution ratio, for pan-sharpening, is typically 4.
%   sensorInf: See generateDefaultSensorInf().
%   Opts: See init_BP_options().
% Output: 
%   I_F_BP: BP optimized sharpening results.
function I_t = BP_Wrapper(I_0, I_MS_LR, I_PAN, ratio, sensorInf, Opts)

is_iter = Opts.is_iter;
is_SS = Opts.is_SS;
DT_mode = Opts.DT_mode;
gamma = Opts.gamma;
mu = Opts.mu;
m_lambda = Opts.m_lambda;
auto_m_lambda = Opts.auto_m_lambda;
if is_iter
    n_reg = Opts.n_reg;
    if is_SS
        RT_mode = Opts.RT_mode;
        I_t = SSBP(I_0, I_MS_LR, I_PAN, ratio, sensorInf, DT_mode, RT_mode, m_lambda, gamma, n_reg);
    else
        if DT_mode
            % BP_T
            interpFunc = @(Er_LR) interpByKernel(Er_LR,ratio,sensorInf.PSF_G.*gamma,sensorInf.upsampling.offset);
        else
            % BP_I
            interpFunc = @(Er_LR) interpWrapper(Er_LR,ratio, sensorInf.upsampling).*gamma;
        end
        I_t = I_0;
        for j = 1:n_reg
            I_t_LR = MTF_conv_sample(I_t, sensorInf, ratio, 1);
            Er_LR = I_MS_LR - I_t_LR;
            Er = interpFunc(Er_LR);
            I_t = I_t + Er;  
        end
        
    end
else
    if is_SS
        RT_mode = Opts.RT_mode;
        Res_mode = 1;
        if auto_m_lambda 
            I_MS = interpWrapper(I_MS_LR,ratio, sensorInf.upsampling);
            I_t1 = FSSBP(I_0, I_PAN, I_MS_LR, ratio, sensorInf, mu, 0, DT_mode, RT_mode, Res_mode, gamma);
            I_t2 = FSSBP(I_0, I_PAN, I_MS_LR, ratio, sensorInf, mu, 1, DT_mode, RT_mode, Res_mode, gamma);
            QP1 = QNR_Plus(I_t1, I_MS_LR, I_MS, I_PAN, 32, sensorInf, ratio);
            QP2 = QNR_Plus(I_t2, I_MS_LR, I_MS, I_PAN, 32, sensorInf, ratio);
            if QP1 > QP2
                lambda = 0;
                I_t = I_t1;
            else
                lambda = 1;
                I_t = I_t2;
            end
            disp([num2str(QP1) '_' num2str(QP2) '_lambda:' num2str(lambda)]);
        else
            I_t = FSSBP(I_0, I_PAN, I_MS_LR, ratio, sensorInf, mu, m_lambda, DT_mode, RT_mode, Res_mode, gamma);
        end
    else
        % FBP
        I_0_LR = MTF_conv_sample(I_0, sensorInf, ratio, 1);
        I_t = BPCloseForm(I_MS_LR-I_0_LR, 0, ratio, sensorInf, mu, DT_mode, gamma)+I_0;
    end
end