function [X,S,Y,e] = t4_svd(Z,rmax,tol)
%T4_SVD   Singular value decomposition of 4-tensor.
%
%    z1 [-------]              [-------] z2       z1 [-------] z2
%   ----|   X   |----( S )-----|   Y   |----  =  ----|   Z   |----
%       [-------]              [-------]             [-------]
%           |                      |                   |   |
%           | n1                   | n2             n1 |   | n2
%
%   See also T4_TT, T4_TMT, T4_MTMTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[z1,z2,n1,n2] = size(Z);

%% svd
Z = permute(Z,[1,3,2,4]);
Z = reshape(Z,z1*n1,z2*n2);
[U,S,V] = svd(Z,'econ');

%% rank
if nargin < 2
    r = rank(S);
elseif nargin < 3
    r = min(rmax,rankdiag(S));
else
    r = min(rmax,rankdiag(S,tol));
end

%% error
if nargout > 3
    if r == size(S,1)
        e = 0;
    else
        e = S(r+1,r+1);
    end
end

%% output
X = reshape(U(:,1:r),z1,n1,r);
X = permute(X,[1,3,2]);
S = S(1:r,1:r);
Y = reshape(V(:,1:r),z2,n2,r);
Y = conj(permute(Y,[3,1,2]));

end


function r = rankdiag(D,tol)

s = diag(D);
if nargin == 1
   tol = max(size(D)) * eps(s(1));
end
r = sum(s > tol);

end
