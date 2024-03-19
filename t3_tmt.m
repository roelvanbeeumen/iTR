function Z = t3_tmt(X,M,Y)
%T3_TMT   3-tensor times matrix times 3-tensor.
%
%   Z = t3_tmt(X,M,Y) returns the 3-tensor Z, which is the tensor product
%   of the 3-tensor X, the matrix M, and the 3-tensor Y.
%
%    x1 [-------] y2       x1 [-----] x2   m1       m2   y1 [-----] y2
%   ----|   Z   |----  =  ----|  X  |---------( M )---------|  Y  |----
%       [-------]             [-----]                       [-----]
%           |                    |                             |
%           | n1*n2              | n1                          | n2
%
%   See also t4_TMT, T3_TT, T3_MTMT, T3_TMTM, T3_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);
[m1,m2] = size(M);     assert(x2 == m1);
[y1,y2,n2] = size(Y);  assert(m2 == y1);

%% Z(:,:,ij) = X(:,:,i)*M*Y(:,:,j)
Z = zeros(x1,y2,n1*n2);
for j = 1:n2
    for i = 1:n1
        idx = i + n1*(j-1);
        Z(:,:,idx) = X(:,:,i)*M*Y(:,:,j);
    end
end

end
