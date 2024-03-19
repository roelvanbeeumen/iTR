function [X,Y] = iTR2(r1,r2,n1,n2)
%ITR2   Infinite tensor ring with 2 cores.
%
%   r1  [-------]  r2     r2  [-------]  r1
%   ----|   X   |---- --- ----|   Y   |----
%       [-------]             [-------]
%           |                     |
%           | n1                  | n2
%
%   See also ITR2C, ITR1, ITR1C.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 2, r2 = r1; end
if nargin < 3, n1 = 2; end
if nargin < 4, n2 = n1; end

%% iTR
X = zeros(r1,r2,n1);
for i = 1:n1
    X(:,:,i) = randn(r1,r2);
end
Y = zeros(r2,r1,n2);
for i = 1:n2
    Y(:,:,i) = randn(r2,r1);
end

end
