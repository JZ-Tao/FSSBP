function [Z] = Sylvester(C1,FBm, FP, ratio,n_row_LR,n_col_LR,C3)
n_band = size(C1,2);
n_row_HR = ratio*n_row_LR;
n_col_HR = ratio*n_col_LR;


% FBm   = fft2(B); 
%FBmC  = FG; %conj(FBm);
if size(FBm,3) == 1
    FBs  = repmat(FBm,[1 1 n_band]);
    FBCs = repmat(FP,[1 1 n_band]);
else
    FBs = FBm;
    FBCs = FP;
end
FBG = FBs.*FBCs;

[Q,Lambda]=eig(C1);
Lambda=reshape(diag(Lambda),[1 1 n_band]);
%InvLbd=1./repmat(Lambda,[ratio*n_row_LR ratio*n_col_LR 1]);
B2Sum=PPlus(FBG./(ratio^2),n_row_LR,n_col_LR);
% B2Sum2=BlockMM(n_row_LR,n_col_LR,ratio^2,n_row_LR*n_col_LR,n_band, FBG);
% B2Sum3 = B2Sum(1:n_row_LR,1:n_col_LR,:);
% Berr = B2Sum2 - B2Sum3;
% norm(Berr(:), 'fro')
%B2Sum=PPlus(abs(FBs).^2./(ratio^2),n_row,n_col);
InvDI=1./(B2Sum(1:n_row_LR,1:n_col_LR,:)+repmat(Lambda,[n_row_LR n_col_LR 1]));
C30=fft2(reshape((Q\C3)',[n_row_HR n_col_HR n_band]))./repmat(Lambda,[ratio*n_row_LR ratio*n_col_LR 1]);

temp  = PPlus_s(C30/(ratio^2).*FBs,n_row_LR,n_col_LR); % The operation: temp=1/d*C5bar*Dv
 

invQUF = C30-repmat(temp.*InvDI,[ratio ratio 1]).*FBCs; % The operation: C5bar- temp*(\lambda_j d Im+\Sum_i=1^d Di^2)^{-1}Dv^H)
VXF    = Q*reshape(invQUF,[n_col_HR*n_col_HR n_band])';
Z = real(ifft2(reshape(VXF',[n_row_HR n_col_HR n_band]))); 


% clip 4 patches, and then sum
% function x = BlockMM(nr,nc,sf,m,n_band, x1)
% myfun = @(block_struct) reshape(block_struct.data,[m,n_band]);
% x1 = blockproc(x1,[nr nc],myfun);
% 
% x1 = reshape(x1,[m,sf,n_band]);
% x1 = mean(x1,2);
% x = reshape(x1,[nr,nc,n_band]);