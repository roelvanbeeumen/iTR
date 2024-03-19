function X = iTR1(r,n)
%ITR1   Infinite tensor ring with 1 core.
%
%    r  [-------]  r
%   ----|   X   |----
%       [-------]
%           |
%           | n
%
%   See also ITR1C, ITR2, ITR2C.

%   Roel Van Beeumen
%   March 18, 2024

%% default argument
if nargin < 2, n = 2; end

%% iTR
X = zeros(r,r,n);
for i = 1:n
    X(:,:,i) = randn(r,r);
end

end
