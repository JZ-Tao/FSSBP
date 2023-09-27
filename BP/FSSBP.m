function I_t = FSSBP(I_0, I_PAN, I_MS_LR, ratio, sensorInf, mu, tau, GT_mode, RT_mode, Res_mode, gamma)
% I_PAN: Z
% I_HS_LR: Y
if ~exist('GT_mode', 'var')
    GT_mode = 1;
end
if ~exist('RT_mode', 'var')
    RT_mode = 1;
end
if ~exist('Res_mode', 'var')
    Res_mode = 0;
end
if ~exist('mu', 'var')
    mu = 1e-3;
end
if ~exist('tau', 'var')
    tau = 1;
end
%max_v = 2^double(L)-1; 
max_v = 1;
I_MS_LR = I_MS_LR./max_v;
I_PAN = I_PAN./max_v;
I_0 = I_0./max_v;

P = I_PAN;
Y = I_MS_LR;

[n_r0, n_c0, ~] = size(P);

P_LR = MTF_conv_sample(P, sensorInf, ratio, 1);
if RT_mode == 1
    [R,~] = estR_new(Y,P_LR);
    P = P - R(1,end);
    P(P<0) = 0;
   
    R = R(:,1:end-1);
    Rt = @ (x) mat2im(R'*im2mat(x), size(x,1));
    RtR = R'*R;
else
    h = estimation_alpha(Y,P_LR,'local');
    R = h';
    
    I_HS = interpWrapper(Y,ratio, sensorInf.upsampling);
    I = mat2im(R*im2mat(I_HS), size(I_HS,1));
    
    for i = 1:size(I_HS,3)
        m = I_HS(:,:,i);
        x = I;
        c = cov(x(:),m(:));
        Gain(i,1) = c(1,2)/var(x(:));
    end    
    Rt = @ (x) mat2im(Gain*im2mat(x), size(x,1));
    RtR = Gain*R;
end


if Res_mode
    P = P - mat2im(R*im2mat(I_0),size(I_0,1));
    Y = Y - MTF_conv_sample(I_0, sensorInf, ratio, 1);
    X_tilde_out = I_0;
    X_tilde_in = 0;
else
    X_tilde_out = 0;
    X_tilde_in = I_0;
end

[n_row, n_col, ~] = size(P);
n_band = size(Y,3);
K_B = sensorInf.PSF_G;

if GT_mode == 1
    K_P = K_B;
else
    gamma = gamma/(ratio^2); % assert: ratio^2 == sum(K_P(:)) == 16
    K_P = getInterpKernel(ratio, sensorInf.upsampling.interp_type, 2, sensorInf.upsampling.tap);
end
    K_P = K_P.*gamma;
usingGGt = 1;
if usingGGt

    G = @(x) MTF_conv_sample(x, sensorInf, ratio, 1);
    if GT_mode
        % assert: sum(K_P(:)) == 1
        Gt = @(x) interpByKernel(x,ratio,K_P,sensorInf.upsampling.offset);
    else
        Gt = @(x) interpWrapper(x, ratio, sensorInf.upsampling).*gamma;
    end

    GGt       = constructGGt_mod(ratio, sensorInf.PSF_G, K_P, n_row,n_col);
    GtY = Gt(Y);
else
    if GT_mode
        FB = zeros(n_row, n_col, size(K_B,3));
        for i = 1:size(K_B,3)
            FB(:,:,i) = psf2otf(K_B(:,:,i), [n_row, n_col]);
        end
        FP = FB;
    else
        FP = psf2otf(K_P, [n_row, n_col]);
    end
    offset = sensorInf.downsampling.offset;
    GtY = zeros(size(I_0));
    GtY(1+offset:ratio:end,1+offset:ratio:end,:) = Y;
    GtY = imfilter(GtY, K_P, 'circular');
end
RtP = Rt(P);
H1 = tau*(RtR) + mu*eye(size(R,2));
H3 = im2mat(tau*RtP + GtY + mu.*X_tilde_in);

if usingGGt
    [Q,Lambda]=eig(H1);
    Lambda=reshape(diag(Lambda),[1 1 n_band]);

    rhs = mat2im((Q\H3),n_row)./repmat(Lambda,[n_row n_col 1]);
    
    DI = GGt + repmat(Lambda,[size(GGt,1) size(GGt,2) 1]);
    I_t = real(mat2im(Q*im2mat(rhs - Gt(ifft2(fft2(G(rhs))./DI))),n_row));
else
    I_t = Sylvester(H1, FB, FP, ratio, size(I_MS_LR,1), size(I_MS_LR,2), H3);
end
I_t = (I_t+X_tilde_out).*max_v;

end