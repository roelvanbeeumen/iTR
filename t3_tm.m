function Y = t3_tm(X,M)
%T3_TM  3-tensor times matrix.
%
%   Y = t3_tm(X,M) returns the 3-tensor Y, which is the tensor product of
%   the 3-tensor X and the matrix M.
%
%    r1 [-----] m2       r1 [-----] r2   m1       m2
%   ----|  Y  |----  =  ----|  X  |---------( M )----
%       [-----]             [-----]
%          |                   |
%          | n                 | n
%
%   See also T3_MT, T3_MTM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[r1,r2,n] = size(X);
[m1,m2] = size(M);
assert(r2 == m1);

%% Y(i) = X(i)*M
Y = zeros(r1,m2,n);
for i = 1:n
    Y(:,:,i) = X(:,:,i)*M;
end

end
