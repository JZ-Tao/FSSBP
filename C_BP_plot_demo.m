% Verification of BP convergence conditions for common sharpening
% interpolators.
% FSSBP: Fast Spatial-Spectral Back Projection Based on Pan-sharpening
% Iterative Optimization, 2023. https://doi.org/10.3390/rs15184543.
sz_X = 56;
% The corresponding convolutional kernel h1 varies somewhat from sensor to
% sensor.
sensor = 'IKONOS'; %'IKONOS' 'GeoEye1' 'WV2' 'WV3','WV3_4bands'  'none'
interplator = 'general';
ratio = 4;

sampling_opts = get_sampling_pars(ratio, interplator);
interp_type = sampling_opts.up.interp_type;
using_default_interp = interp_type ~= 1;
tap = sampling_opts.up.tap;
N = 41;
fcut = 1/ratio;
PSF_G = zeros(N,N,1);
GNyqs = 0.1:0.05:0.5;
for jj = 1:length(GNyqs)
    GNyq = GNyqs(jj);
    for ii = 1 : 1
        alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq(ii))));
        H = fspecial('gaussian', N, alpha);
        Hd = H./max(H(:));
        h = fwind1(Hd,kaiser(N));
        PSF_G(:,:,ii) = real(h);
    end
    h1 = PSF_G./sum(PSF_G(:));

    % Convolution kernel h2 (no need to focus on details). 
    % The shape is similar to a Gaussian surface, but non-Gaussian.
    if using_default_interp 
        Kernel_scaling = 2;
    else
        Kernel_scaling = 4;
    end
    h2 = getInterpKernel(Kernel_scaling, interp_type, 2, tap);
    h3 = getInterpKernel(Kernel_scaling, interp_type, 2, 7);
    h4 = getInterpKernel(Kernel_scaling, interp_type, 2, 8);
    h5 = getInterpKernel(Kernel_scaling, interp_type, 2, 23);
    X = randn(sz_X,sz_X);%ones(8,8);


    x = X(:);
    h = h1;
    p = h2;
    [ggt_f,ngs_t(jj)] = constructHP(h,h,ratio,sz_X,sz_X);
    [ggt_f2,ngs_gen(jj)] = constructHP(h,h2,ratio,sz_X,sz_X);
    [ggt_f3,ngs_tap7(jj)] = constructHP(h,h3,ratio,sz_X,sz_X);
    [ggt_f5,ngs_tap23(jj)] = constructHP(h,h5,ratio,sz_X,sz_X);
end


figure, plot(GNyqs,ngs_t,'--',GNyqs,ngs_gen,':',GNyqs,ngs_tap7,'.-',GNyqs,ngs_tap23,'-'); title('Verification of convergence condition');
legend('MTF','general','tap-7','tap-23');
xlabel('MTF Nyquist Frequency');
ylabel('C_{BP}');
