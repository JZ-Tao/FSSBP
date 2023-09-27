function P = generatePaddingMatrix(sz_x, sz_p)
% Equivalent to a boundary extension of X.
% Input image size sz_x, and extension size sz_p.
% Idea: In P, whenever there is an element in the target PX matrix that is
% not 0, there is a 1 in the row corresponding to the index number of that
% element.
% And for each such row, the value of 1 is shifted 1 column to the right.

% TODO: Loop Boundary Handling
PX = padarray(ones(sz_x), sz_p);
indx = find(PX~=0);
px = (sz_x(1)+sz_p(1)*2)*(sz_x(2)+sz_p(2)*2);
py = sz_x(1)*sz_x(2);
P = zeros(px,py);
for i=1:length(indx)
    P(indx(i),i) = 1;
end


