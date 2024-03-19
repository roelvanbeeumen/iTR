function Z = t3_mtmtm(A,X,B,Y,C)
%T3_MTMTM   Matrix times 3-tensor times matrix times 3-tensor times matrix.
%
%   Z = t3_mtmtm(A,X,B,Y,C) returns the 3-tensor Z, which is the tensor
%   product of the matrix A, the 3-tensor X, the matrix B, the 3-tensor Y,
%   and the matrix C.
%
%    a1 [-------] c2       a1       a2   x1 [-----] x2   b1       b2   y1 [-----] y2   c1       c2
%   ----|   Z   |----  =  ----( A )---------|  X  |---------( B )---------|  Y  |---------( C )----
%       [-------]                           [-----]                       [-----]
%           |                                  |                             |
%           | n1*n2                            | n1                          | n2
%
%   See also T4_MTMTM, T3_TT, T3_TMT, T3_MTMT, T3_TMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[a1,a2] = size(A);
[x1,x2,n1] = size(X);  assert(a2 == x1);
[b1,b2] = size(B);     assert(x2 == b1);
[y1,y2,n2] = size(Y);  assert(b2 == y1);
[c1,c2] = size(C);     assert(y2 == c1);

%% Z(:,:,ij) = A*X(:,:,i)*B*Y(:,:,j)*C
Z = zeros(a1,c2,n1*n2);
for j = 1:n2
    for i = 1:n1
        idx = i + n1*(j-1);
        Z(:,:,idx) = A*X(:,:,i)*B*Y(:,:,j)*C;
    end
end

end
