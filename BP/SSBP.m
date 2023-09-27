function I_t = SSBP(BP_I0, I_MS_LR, I_PAN, ratio, sensorInf, DT_mode, RT_mode, m_lambda, gamma, n_reg)
I_t = BP_I0;
if isempty(BP_I0) 
    disp('Empty BP_init. Skipping SSBP');
    return;
end
if RT_mode == 1
    [R,~] = estR_new(I_MS_LR, MTF_conv_sample(I_PAN, sensorInf, ratio, 1));
    P = I_PAN - R(1,end);
    P(P<0) = 0;
    R = R(:,1:end-1);
    F_Rt = @ (x) mat2im(R'*im2mat(x), size(x,1));
    F_R = @(x) mat2im(R*im2mat(x), size(x,1));
elseif RT_mode == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    P = I_PAN;
    P_LR = MTF_conv_sample(P, sensorInf, ratio,1);
    I_MS = interpWrapper(I_MS_LR,ratio, sensorInf.upsampling);
    P_LP = interpWrapper(P_LR,ratio, sensorInf.upsampling);
    h = estimation_alpha(I_MS_LR,P_LR,'local');
    alpha(1,1,:) = h;
    F_R = @(x) sum(x.* repmat(alpha,[size(x,1) size(x,2) 1]),3);

    Gain = zeros(size(I_MS,3),1);
    for gi = 1:size(I_MS,3)
        m = I_MS(:,:,gi);
        x = P_LP;%I;
        c = cov(x(:),m(:));
        Gain(gi,1) = c(1,2)/var(x(:));
    end    
    F_Rt = @ (x) mat2im(Gain*im2mat(x), size(x,1));%Gain.*repmat(x,[1,1,size(Y,3)]);
end

if DT_mode
    interpFunc = @(Er_LR) interpByKernel(Er_LR,ratio,sensorInf.PSF_G.*gamma,sensorInf.upsampling.offset);
else
    interpFunc = @(Er_LR) interpWrapper(Er_LR,ratio, sensorInf.upsampling).*gamma;
end
%I_t = I_0;
%m_lambda = 0;
for j = 1:n_reg
    I_t_LR = MTF_conv_sample(I_t, sensorInf, ratio, 1);
    Er_LR = I_MS_LR - I_t_LR;
    R_S = interpFunc(Er_LR);
    R_lambda = F_Rt(P-F_R(I_t));
    R_SS = R_S + m_lambda.*R_lambda;
    I_t = I_t + R_SS;
end
