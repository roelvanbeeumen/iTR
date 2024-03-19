function Z = t4_tt(X,Y)
%T4_TT  3-tensor times 3-tensor.
%
%   Z = t4_tt(X,Y) returns the 4-tensor Z, which is the tensor product of
%   the 3-tensors X and Y.
%
%    x1 [-------] y2       x1 [-----] x2   y1 [-----] y2
%   ----|   Z   |----  =  ----|  X  |---------|  Y  |----
%       [-------]             [-----]         [-----]
%         |   |                  |               |
%      n1 |   | n2               | n1            | n2
%
%   See also T3_TT, T4_TMT, T4_MTMT, T4_TMTM, T4_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);
[y1,y2,n2] = size(Y);
assert(x2 == y1);

%% Z(:,:,i,j) = X(:,:,i)*Y(:,:,j)
Z = zeros(x1,y2,n1,n2);
for i = 1:n1
    for j = 1:n2
        Z(:,:,i,j) = X(:,:,i)*Y(:,:,j);
    end
end

end
