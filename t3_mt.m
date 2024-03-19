function Y = t3_mt(M,X)
%T3_MT  Matrix times 3-tensor.
%
%   Y = t3_mt(M,X) returns the 3-tensor Y, which is the tensor product of
%   the matrix M and the 3-tensor X.
%
%    m1 [-----] r2       m1       m2   r1 [-----] r2
%   ----|  Y  |----  =  ----( M )---------|  X  |----
%       [-----]                           [-----]
%          |                                 |
%          | n                               | n
%
%   See also T3_TM, T3_MTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[m1,m2] = size(M);
[r1,r2,n] = size(X);
assert(m2 == r1);

%% Y(i) = M(i)*X
Y = zeros(m1,r2,n);
for i = 1:n
    Y(:,:,i) = M*X(:,:,i);
end

end
