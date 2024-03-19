function Z = t3_tt(X,Y)
%T3_TT  3-tensor times 3-tensor.
%
%   Z = t3_tt(X,Y) returns the 3-tensor Z, which is the tensor product of
%   the 3-tensors X and Y.
%
%    x1 [-------] y2       x1 [-----] x2   y1 [-----] y2
%   ----|   Z   |----  =  ----|  X  |---------|  Y  |----
%       [-------]             [-----]         [-----]
%           |                    |               |
%           | n1*n2              | n1            | n2
%
%   See also T4_TT, T3_TMT, T3_MTMT, T3_TMTM, T3_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);
[y1,y2,n2] = size(Y);
assert(x2 == y1);

%% Z(:,:,ij) = X(:,:,i)*Y(:,:,j)
Z = zeros(x1,y2,n1*n2);
for j = 1:n2
    for i = 1:n1
        idx = i + n1*(j-1);
        Z(:,:,idx) = X(:,:,i)*Y(:,:,j);
    end
end

end
