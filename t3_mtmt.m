function Z = t3_mtmt(A,X,B,Y)
%T3_MTMT  Matrix times 3-tensor times matrix times 3-tensor.
%
%   Z = t3_mtmt(A,X,B,Y) returns the 3-tensor Z, which is the tensor
%   product of the matrix A, the 3-tensor X, the matrix B, and the
%   3-tensor Y.
%
%    a1 [-------] y2       a1       a2   x1 [-----] x2   b1       b2   y1 [-----] y2
%   ----|   Z   |----  =  ----( A )---------|  X  |---------( B )---------|  Y  |----
%       [-------]                           [-----]                       [-----]
%           |                                  |                             |
%           | n1*n2                            | n1                          | n2
%
%   See also T4_MTMT, T3_TT, T3_TMT, T3_TMTM, T3_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[a1,a2] = size(A);
[x1,x2,n1] = size(X);  assert(a2 == x1);
[b1,b2] = size(B);     assert(x2 == b1);
[y1,y2,n2] = size(Y);  assert(b2 == y1);

%% Z(:,:,ij) = A*X(:,:,i)*B*Y(:,:,j)
Z = zeros(a1,y2,n1*n2);
for j = 1:n2
    for i = 1:n1
        idx = i + n1*(j-1);
        Z(:,:,idx) = A*X(:,:,i)*B*Y(:,:,j);
    end
end

end
