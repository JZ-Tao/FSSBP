%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: 
%           interp23tap interpolates the image I_Interpolated using a polynomial with 23 coefficients interpolator. 
% 
% Interface:
%           I_Interpolated = interp23tap(I_Interpolated,ratio)
%
% Inputs:
%           I_Interpolated: Image to interpolate;
%           ratio:          Scale ratio between MS and PAN. Pre-condition: Resize factors power of 2;
%
% Outputs:
%           I_Interpolated: Interpolated image.
% 
% References:
%           [Aiazzi02]      B. Aiazzi, L. Alparone, S. Baronti, and A. Garzelli, Context-driven fusion of high spatial and spectral resolution images based on
%                           oversampled multiresolution analysis,?IEEE Transactions on Geoscience and Remote Sensing, vol. 40, no. 10, pp. 2300?312, October
%                           2002.
%           [Aiazzi13]      B. Aiazzi, S. Baronti, M. Selva, and L. Alparone, Bi-cubic interpolation for shift-free pan-sharpening,?ISPRS Journal of Photogrammetry
%                           and Remote Sensing, vol. 86, no. 6, pp. 65?6, December 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function I_Interpolated = interpDefault(I_Interpolated,ratio,tap, offset)
if ~exist('tap','var')
    tap = -1;
end
if ~exist('offset','var')
    offset = [1,0];
end
if (2^round(log2(double(ratio))) ~= ratio)
    disp('Error: Only resize factors power of 2');
    return;
end 

[r,c,b] = size(I_Interpolated);
BaseCoeff = getInterpKernel(2, 2, 1, tap);

first = 1;
%BaseCoeff2 = conv2(BaseCoeff',BaseCoeff);
% ksz = floor(size(BaseCoeff,2)/2);
for z = 1 : log2(ratio)

    I1LRU = zeros((2^z) * r, (2^z) * c, b);
    
    if first
        I1LRU(1+offset(1):2:end,1+offset(1):2:end,:) = I_Interpolated;
        %I1LRU(2:2:end,2:2:end,:) = I_Interpolated;
        first = 0;
    else
        I1LRU(1+offset(2):2:end,1+offset(2):2:end,:) = I_Interpolated;
    end

     for ii = 1 : b
        t = I1LRU(:,:,ii); 
        %I1LRU2(:,:,ii) = imfilter(t,BaseCoeff2,'circular'); 
        t = imfilter(t',BaseCoeff,'circular'); 
        I1LRU(:,:,ii) = imfilter(t',BaseCoeff,'circular'); 
%         t = imfilter(t',BaseCoeff); 
%         I1LRU(:,:,ii) = imfilter(t',BaseCoeff); 
    end

    I_Interpolated = I1LRU;
    
end

end