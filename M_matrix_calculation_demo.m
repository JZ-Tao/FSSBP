% Demonstration of the BP transformation matrix and its spectral radius
% calculation.
% FSSBP: Fast Spatial-Spectral Back Projection Based on Pan-sharpening
% Iterative Optimization, 2023. https://doi.org/10.3390/rs15184543.
% =========================================================================
% Note: The theoretical problem is to apply a transformation M to a
% two-dimensional matrix X to obtain an output y, i.e., y = Mx.
% where x is the vectorized representation of X, and M is the matrix form
% corresponding to a series of "operations" (e.g., convolution, upsampling,
% downsampling).
% The operations applied are as follows (in order from left to right):
% Convolution h1, 4x downsampling, 2x upsampling, convolution h2, 2x
% upsampling, convolution h2.
% Theoretically corresponding matrix M = H2*U2*H2*U1*S*H1 (left
% multiplication by x).
% =========================================================================
% The following program flow consists of two parts.
% The first part consists only of the "operations" and does not give M
% explicitly.
% The second part illustrates the explicit exploration of the contents of
% the matrix of M.
% =========================================================================

% Part 1: Example of Y = f(X) operation.
% Program input:
% X: a two-dimensional matrix of arbitrary size in multiples of 4 (default
% 16*16, should not exceed 64*64) with arbitrary content.
% h1: convolution kernel 1 (size 41 * 41, Gaussian low-pass filtering), the
% corresponding matrix can be written as H1.
% h2: convolution kernel 2 (size 23*23, low-pass filtering), the
% corresponding matrix can be written as H2.
% Program output:
% Y: processed 2D matrix, size equal to X.
% =========================================================================

% =========================================================================
% Second part of the program: Generate M based on x dimension information
% and given convolution kernel and sampling information.
% 
% Program input:
% x: vectorization of X
% h1: convolution kernel 1 (size 41*41, Gaussian low-pass filtering), the
% corresponding matrix form can be written as H1.
% h2: convolution kernel 2 (size 23*23, low-pass filtered), the
% corresponding matrix form can be written as H2.
% Program output:
% y2: y2 = M*x. ideally with y2 == Y(:).
% =========================================================================
% Author: jingzhe tao. 2020.7
% =========================================================================
close all; clear;
%% Input data preparation
% Adjustable parameters
sz_X = 56;
enablePart2 = 1;
enablePart3 = 0;
enablePart4 = 0;
% The convolution kernel h1 corresponding to different sensors varies
% somewhat.
sensor = 'IKONOS'; %'IKONOS' 'GeoEye1' 'WV2' 'WV3','WV3_4bands'  'none'
interplator = 'general';
ratio = 4;
   
% Convolution kernel h1. The convolution kernel is a 2D Gaussian surface
% shape.
sampling_opts = get_sampling_pars(ratio, interplator);
tap = sampling_opts.up.tap;
sensorInf.upsampling = sampling_opts.up;
sensorInf.downsampling = sampling_opts.down;
% If the first offset is non-zero, after a second 2-fold upsampling, an
% increment will be generated, e.g., [1,0] will correspond to offset==2
% (i.e., corresponds to subscript == 3).
offset_half = sampling_opts.up.offset;
offset = sampling_opts.down.offset;
% using_imresize = 0;
using_default_interp = sensorInf.upsampling.interp_type ~= 1;

h1 = Blur_Kernel(1, sensor, ratio);

N = 17;
fcut = 1/ratio;
PSF_G = zeros(N,N,1);
GNyq = 0.2;
for ii = 1 : 1
    alpha = sqrt((N*(fcut/2))^2/(-2*log(GNyq(ii))));
    H = fspecial('gaussian', N, alpha);
    Hd = H./max(H(:));
    h = fwind1(Hd,kaiser(N));
    PSF_G(:,:,ii) = real(h);
end
h1 = PSF_G./sum(PSF_G(:));

% Convolution kernel h2. The shape is similar to a Gaussian surface, but
% non-Gaussian.
if using_default_interp 
    Kernel_scaling = 2;
else
    Kernel_scaling = 4;
end
h2 = getInterpKernel(Kernel_scaling, sensorInf.upsampling.interp_type, 2, tap);
%h2 = h1; %%%%%%%%%%%%%%%%%%%%

%rng(2333);
X = randn(sz_X,sz_X);%ones(8,8);

% X = double(rgb2gray((imread('image/GT3N.png'))));
% X = padarray(X,[100,100]);
% sz_X = size(X,1);

% X = otf2psf(psf2otf([0 0 0; 0 1 0; 0 0 0], [sz_X,sz_X]));
if sz_X >= 60
    warning('X size is too large leading to slow computation.'); 
    warning('Explicit construction of matrix M will be turned off, other parts are not affected.');
    enablePart2 = 0;
end
% if sz_X < size(h1,1)
%     warning('X size is smaller than convolutional kernel.');
%     warning('OTF conversion does not work and is temporarily disabled. part3 and 4 will be closed, other parts are not affected');
%     enablePart3 = 0;
%     enablePart4 = 0;
% end
x = X(:);
h = h1;
p = h2;
ggt_f = constructHP(h,p,ratio,sz_X,sz_X);
% h_p = ifft2(ggt_f);
% norm(h_p, 1)

% h_p = conv2(h, p);
% h_p = imresize(h_p, ratio, 'nearest');
%h_p = downsampleWrapper(h_p, ratio, using_imresize);
%mid_hp = round((size(h_p,1)+1)/2);
%h_p(mid_hp, mid_hp) = 1- h_p(mid_hp, mid_hp);

%norm(h_p, 2)
% delta = zeros(size(h_p));
% delta(mid_hp, mid_hp) = 1;
% delta_h = conv2(delta, h);
% h_p2 = conv2(delta_h, p);
%% Part I
disp('Part I: Schematic execution of an image processing operation');
% Operation 1: image filtering (convolution), theoretically with
% corresponding matrix multiplication form H1x = H1*x.
X1 = imfilter(X,h1,'circular');
% Operation 2: 4-fold downsampling, i.e., take 1 every 4 for X1 elements.
% SH1x = S*H1x.
X2 = downsampleWrapper(X1, ratio, sensorInf.downsampling);
%X2 = downsample(downsample(X1,ratio,1)',ratio,1)';%X2 = imresize(X1, 1/ratio, 'nearest');
if using_default_interp
    % Operation 3.1A: 2x upsampling, i.e., adding 0, to X2.
    X3 = zeros(ratio/2*size(X2,1), ratio/2*size(X2,2));    
    X3(1+offset_half(1):2:end,1+offset_half(1):2:end,:) = X2;
    % Operation 3.2A: Image Filtering (Convolution).
    X4 = imfilter(X3,h2); %X4 = imfilter(X3,h2,'circular');
    % Operation 3.3A: 2x upsampling, i.e., adding 0 to X2
    X5 = zeros(ratio*size(X2,1), ratio*size(X2,2));
    X5(1+offset_half(2):2:end,1+offset_half(2):2:end,:) = X4;
else
    % Operation 3B: 4x upsampling.
    X5 = zeros(ratio.*size(X2,1), ratio.*size(X2,2));
    X5(1+offset:ratio:end,1+offset:ratio:end) = X2;
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X5 = zeros(size(X1));
X5(1+offset:ratio:end,1+offset:ratio:end) = X1(1+offset:ratio:end,1+offset:ratio:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operation 4: Image Filtering (Convolution).
X6 = imfilter(X5,h2,'circular'); 

% Output Y.
Y = X6;
% Y = interpWrapper(X2, ratio, sensorInf.upsampling);
% = ===================================================================== =
% Note that the goal is not to derive Y, but rather the transformation
% matrix corresponding from X --> Y, denoted M.
% implies that Y = f(X) <----> y = Mx.
% In the theoretical derivation, M will be relatively abstract and the
% details of its matrix values may not be known. Different x's correspond to different dimensions of M.
% = ===================================================================== =

%% Part II
if enablePart2
    % The following matrix equivalence forms are explored for each
    % operation.
    disp('Part II: Explicit Construction of Matrix M');
    % Step 1
    % H1, note that the reason the matrix size is not square is because
    % there is a bounding extension (in this case 0's complement, the
    % actual circular convolution used is not 0's complement).
    % The boundary extension can also be theoretically written in the form
    % of matrix multiplication, i.e., (H1P)x, with P being the extension matrix.
    % The extension is: A = [1 2; 3 4], B = padarray(A, [1 1]), a = A(:), b
    % = B(:), then Pa = b, where P is of size 16 * 4.
    % Let X be of size m*n, x be u*1 (u=m*n),h1 be v*v, then the size of H1
    % be u*[(m+v-1)*(n+v-1)].
    % The general principle is that the convolution operation can be
    % transformed into a Toplitz matrix, but for practical boundary
    % extension processing reasons, resulting in H1 dimensions varying from
    % x to x, roughly as a chunked matrix repeated along the diagonal
    % direction (banded matrices?).
    [~, H1] = conv2tp2d(X, h1);
    % Here the extension matrix P is taken into account, implicitly in H1,
    % so that HP1 = H1*P1.
    P1 = generatePaddingMatrix(size(X), (size(h1)-1)./2);
    HP = H1*P1;
    %figure, imshow(HP, []); title('The matrix H1 corresponding to the convolution kernel h1');
    % As validation, there should be HP*x == X1(:) <-- imfilter(X,h1);
    T = HP;
    %[X_T, Y_T] = eig(T);
    % disp(['Matrix HP dimensions: ' num2str(size(T,1)) '*'...
    %     num2str(size(T,2)) '; range: ' num2str(min(T(:))) ' - ' num2str(max(T(:)))]);

    % Step 2
    % = ===================================================================== =
    % Further, the next 4-fold sampling after h1 (rows and columns are
    % separated by 4 to take 1), i.e., SH.
    SH = zeros(size(HP,1)/(ratio^2), size(HP,2));
    for i = 1:size(HP,2)
        %a = downsample(downsample(reshape(HP(:,i), size(X)),ratio,1)',ratio,1)'; 
        a = downsampleWrapper(reshape(HP(:,i), size(X)), ratio, sensorInf.downsampling);
        SH(:,i) = a(:);
    end
    HTSTSH = SH'*SH; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    GTG = HTSTSH;
    GGT = SH*SH';
    [~, GTG_lambda] = eig(GTG);
    [~, GGT_lambda] = eig(GGT);
    % As a verification, there should be SH*x == X2(:).
    %figure, imshow(SH, []); title('left multiplied by S, i.e. SH');
    if using_default_interp
        % Step 3-1A
        % Ditto for 2x sampling.
        U1SH = zeros(size(SH,1)*(ratio/2)^2, size(SH,2)); % (ratio/2)^2
        for i = 1:size(SH,2) 
            a = upsample(upsample(reshape(SH(:,i), size(X)./ratio),ratio/2, offset_half(1))',ratio/2, offset_half(1))'; 
            U1SH(:,i) = a(:);
        end

        % Step 3-2A
        % Interpolation, Low-Pass Convolution
        %U1SHx = U1SH*x; % U1SHx == X3(:)
        [~, H2] = conv2tp2d(X3, h2);
        % The extension matrix P is taken into account such that HP2 =
        % H2*P2.
        P2 = generatePaddingMatrix(size(X3), (size(h2)-1)./2);
        HP2 = H2*P2;
        H2USH = HP2*U1SH;
        % Ditto for 2x sampling.
        % Step 3-3A
        U2H2U1SH = zeros(size(H2USH,1)*(ratio/2)^2, size(H2USH,2));
        for i = 1:size(SH,2) 
            a = upsample(upsample(reshape(H2USH(:,i), size(X)./(ratio/2)), ratio/2, offset_half(2))', ratio/2, offset_half(2))'; 
            U2H2U1SH(:,i) = a(:);
        end

    else
        % Step 3-1B
        U2H2U1SH = zeros(size(SH,1)*(ratio^2), size(SH,2));
        for i = 1:size(SH,2)
    %         t = reshape(SH(:,i), size(X)./ratio);
    %         T = zeros(ratio.*size(t));
    %         T(floor(ratio/2)+1:ratio:end,floor(ratio/2)+1:ratio:end) = t;
    %         U1SH_1(:,i) = T(:);
            a = upsample(upsample(reshape(SH(:,i), size(X)./ratio),ratio, offset)',ratio, offset)'; 
            U2H2U1SH(:,i) = a(:);
        end
    end
    % Step 4
    %X5(:) == U2H2U1SHx = U2H2U1SH*x;
    [~, H4] = conv2tp2d(X5, h2);
    % The extension matrix P is taken into account such that HP2 = H2*P2.
    P = generatePaddingMatrix(size(X5), (size(h2)-1)/2);
    HP4 = H4*P;
    M = HP4*U2H2U1SH;

    %% Display, statistics of M
    figure, imshow(M, []); title('The matrix M');

    Y2 = reshape(M*x, size(X));
    disp('|| Y2 - Y ||2:'); norm(Y2 - Y, 'fro')

    [~, M_lambda] = eig(M);
    disp(['The spectral radius of M:' num2str(max(abs(diag(M_lambda))))]); 
    disp(['The trace of M:' num2str(trace(M))]);
    I = eye(size(M));
    W = M+I;
    I_M = I-M;
    [~, I_M_lambda] = eig(I_M);
    disp(['The spectral radius of I-M:' num2str(max(abs(diag(I_M_lambda))))]); 
    disp(['The L1-norm of M:' num2str(norm(M,1))]);
    disp(['The L1-norm of I-M:' num2str(norm(I_M,1))]);
    % Note: When i is large enough, I_MM == IM*IMM ~= 0
%         I_ = eye(sz_X);
%         for i = 1:size(I_,2) 
%             a = upsample(upsample(reshape(I_(:,i), size(X)./ratio),ratio, offset)',ratio, offset)'; 
%             ST(:,i) = a(:);
%         end
%         PST = HP4*ST;
%         SHPST = SH*PST;
%         I_W = eye(size(SHPST)) - SHPST;
%         norm(I_W)
%     [~, I_W_lambda] = eig(I_W);
%     disp(['The L1-norm of I_W' num2str(norm(I_W,1))]);
%     disp(['The spectral radius of I_W:' num2str(max(abs(diag(I_W_lambda))))]);

end


