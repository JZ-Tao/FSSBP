function x = BPCloseForm(y, xtilde, ratio, sensorInf, rho, GT_mode, gamma)
if ~exist('gamma','var')
    gamma = 1;
end
%max_v = max(2^L-1, 1);
max_v = 1;
y = y./max_v;
xtilde = xtilde./max_v;

[n_r0, n_c0, ~] = size(y);

G = @(x) MTF_conv_sample(x, sensorInf, ratio, 1);
if ~GT_mode
    gamma = gamma/(ratio^2);
    K_P = getInterpKernel(ratio, sensorInf.upsampling.interp_type, 2, sensorInf.upsampling.tap).*gamma;
    Gt = @(x) interpWrapper(x, ratio, sensorInf.upsampling).*gamma;
else
    K_P = sensorInf.PSF_G.*gamma;
    Gt = @(x) interpByKernel(x,ratio,K_P,sensorInf.upsampling.offset);
end

Gty       = Gt(y);

y = imresize(Gty,1/ratio);
[n_r,n_c, n_b] = size(y);
rows      = n_r*ratio;
cols      = n_c*ratio;

GGt       = constructGGt_mod(ratio, sensorInf.PSF_G, K_P, rows,cols);

rhs = Gty + rho*xtilde;
x = (rhs - Gt(ifft2(fft2(G(rhs))./(GGt + rho))))/rho;
x = real(x).*max_v;

end
