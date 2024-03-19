function Z = t4_mtmt(A,X,B,Y)
%T4_MTMT  Matrix times 3-tensor times matrix times 3-tensor.
%
%   Z = t4_mtmt(A,X,B,Y) returns the 4-tensor Z, which is the tensor
%   product of the matrix A, the 3-tensor X, the matrix B, and the
%   3-tensor Y.
%
%    a1 [-------] y2       a1       a2   x1 [-----] x2   b1       b2   y1 [-----] y2
%   ----|   Z   |----  =  ----( A )---------|  X  |---------( B )---------|  Y  |----
%       [-------]                           [-----]                       [-----]
%         |   |                                |                             |
%      n1 |   | n2                             | n1                          | n2
%
%   See also T3_MTMT, T4_TT, T4_TMT, T4_TMTM, T4_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[a1,a2] = size(A);
[x1,x2,n1] = size(X);  assert(a2 == x1);
[b1,b2] = size(B);     assert(x2 == b1);
[y1,y2,n2] = size(Y);  assert(b2 == y1);

%% Z(:,:,i,j) = A*X(:,:,i)*B*Y(:,:,j)
Z = zeros(a1,y2,n1,n2);
for i = 1:n1
    for j = 1:n2
        Z(:,:,i,j) = A*X(:,:,i)*B*Y(:,:,j);
    end
end

end
