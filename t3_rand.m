function X = t3_rand(r1,r2,n)
%T3_RAND  Uniformly distributed pseudorandom 3-tensor.
%
%    r1 [-------] r2
%   ----|   X   |----
%       [-------]
%           |
%           | n
%
%   See also T3_RANDN, RAND.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 2, r2 = r1; end
if nargin < 3, n = 2; end

%% iTR
X = zeros(r1,r2,n);
for i = 1:n
    X(:,:,i) = rand(r1,r2);
end

end
