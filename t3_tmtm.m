function Z = t3_tmtm(X,B,Y,C)
%T3_TMTM  3-tensor times matrix times 3-tensor times matrix.
%
%   Z = t3_tmtm(X,B,Y,C) returns the 3-tensor Z, which is the tensor
%   product of the 3-tensor X, the matrix B, the 3-tensor Y, and the
%   matrix C.
%
%    x1 [-------] c2       x1 [-----] x2   b1       b2   y1 [-----] y2   c1       c2
%   ----|   Z   |----  =  ----|  X  |---------( B )---------|  Y  |---------( C )----
%       [-------]             [-----]                       [-----]
%           |                    |                             |
%           | n1*n2              | n1                          | n2
%
%   See also T4_TMTM, T3_TT, T3_TMT, T3_MTMT, T3_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);
[b1,b2] = size(B);     assert(x2 == b1);
[y1,y2,n2] = size(Y);  assert(b2 == y1);
[c1,c2] = size(C);     assert(y2 == c1);

%% Z(:,:,ij) = X(:,:,i)*B*Y(:,:,j)*C
Z = zeros(x1,c2,n1*n2);
for j = 1:n2
    for i = 1:n1
        idx = i + n1*(j-1);
        Z(:,:,idx) = X(:,:,i)*B*Y(:,:,j)*C;
    end
end

end
