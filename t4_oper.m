function Z = t4_oper(H,X,Y)
%T4_OPER  Apply operator to 4-tensor.
%
%   Z = t4_oper(H,X) returns the 4-tensor
%
%        x1 [-------] x2       x1 [-------] x2
%       ----|   Z   |----  =  ----|   X   |----
%           [-------]             [-------]
%             |   |                 |   |
%          n1 |   | n2           n1  \ /  n2
%                                     |
%                                   /---\
%                                  (  H  )
%                                   \---/
%                                     |
%                                    / \
%                                n1 |   | n2
%
%   Z = t4_oper(H,X,Y) returns the 4-tensor
%
%        x1 [-------] y2       x1 [-----] x2   y1 [-----] y2
%       ----|   Z   |----  =  ----|  X  |---------|  Y  |----
%           [-------]             [-----]         [-----]
%             |   |                  |               |
%          n1 |   | n2            n1 |               | n2
%                                    |               |
%                                  /-------------------\
%                                 (          H          )
%                                  \-------------------/
%                                    |               |
%                                 n1 |               | n2
%
%   See also T4_TT, T4_TMT, T4_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[h1,h2] = size(H);
if nargin < 3
    [x1,x2,n1,n2] = size(X);
else
    [x1,x2,n1] = size(X);
    [y1,y2,n2] = size(Y);  assert(x2 == y1);
end
assert(n1*n2 == h1);
assert(n1*n2 == h2);

%% apply operator
if nargin < 3
    Z = reshape(X,x1*x2,n1*n2)*H;
else
    Z = t3_tt(X,Y);
    Z = reshape(Z,x1*y2,n1*n2)*H;
end

%% convert to 4-tensor
Z = reshape(Z,x1,x2,n1,n2);

end
