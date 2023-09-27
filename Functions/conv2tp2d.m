% By TomHeaven, hanlin_tan@nudt.edu.cn @ 2016.09.28
function [z, H] = conv2tp2d(X, Y)

sz_y2 = (size(Y,1)-1)/2;
XO = X; X = padarray(X,[sz_y2 sz_y2]); XT = X'; x = XT(:);

% construct convolutional multiplication matrix H
Y_hat = [ fliplr(flipud(Y)) zeros(size(Y, 1), size(X,2) - size(Y, 2))];
y_hat = reshape(Y_hat', [1, size(Y_hat,1) * size(Y_hat,2)]);
y_hat = [ y_hat zeros(1, (size(X,1) - size(Y,1)) * size(X,2) ) ];
 
H = zeros(size(XO,1) * size(XO,2), length(x));
%H(1, :) = y_hat;
len = length(y_hat);
cnt = 0;
for i = 1: size(XO,1)
    for j = 1 : size(XO,2)
       cnt = cnt + 1;
       H(cnt,:) = y_hat;
       y_hat(2:len) = y_hat(1:len - 1);
       y_hat(1) = 0;
    end
    % skip invalid convolution
    for j = 1 :   size(Y, 2) - 1
       y_hat(2:len) = y_hat(1:len - 1);
       y_hat(1) = 0;
    end
end

 
z = H * x;
z = reshape(z', [(size(X,2) - size(Y,2) + 1), (size(X,1) - size(Y,1) + 1)   ])';
 
Z = conv2(X, Y, 'valid');
ZO = conv2(XO, Y, 'same');
% diff = norm(z - Z)
% diff2 = norm(z - ZO)