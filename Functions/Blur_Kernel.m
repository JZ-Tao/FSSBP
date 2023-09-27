function [PSF_G, GNyq] = Blur_Kernel(n_band, sensor, ratio, is_XS)

if ~exist('is_XS','var')
    is_XS = true;
end

GNyq = getGNyqBySensor(sensor, n_band, is_XS);
tap = 41;
PSF_G = MTF_GNyq2PSF(GNyq, tap, ratio);
