function [GGt, ngs] = constructHP(h,p,K,rows,cols)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eigen-decomposition for super-resolution
% Stanley Chan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_p = conv2(h,p);

yc = ceil(size(h_p,1)/2);  % mark the center coordinate
xc = ceil(size(h_p,2)/2);

L = floor(size(h_p,1)/K);  % width of the new filter 
                           %  = (1/k) with of the original filter               
g = zeros(L,L);            % initialize new filter
for i=-floor(L/2):floor(L/2)
    for j=-floor(L/2):floor(L/2)
        g(i+floor(L/2)+1,j+floor(L/2)+1) = h_p(yc+K*i, xc+K*j);
    end
end
g_yc = ceil(size(g,1)/2);  % mark the center coordinate
g_xc = ceil(size(g,2)/2);
gs(g_yc, g_xc) = 1 - g(g_yc, g_xc);
ngs = norm(gs,1);
GGt = abs(fft2(g,rows/K,cols/K));