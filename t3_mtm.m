function Y = t3_mtm(A,X,B)
%T3_MTM   Matrix times 3-tensor times matrix.
%
%   Y = t3_mtm(A,X,B) returns the 3-tensor Y, which is the tensor product
%   of the matrix A, the 3-tensor X, and the matrix B.
%
%    a1 [-----] b2       a1       a2   r1 [-----] r2   b1       b2
%   ----|  Y  |----  =  ----( A )---------|  X  |---------( B )----
%       [-----]                           [-----]
%          |                                 |
%          | n                               | n
%
%   See also T3_MT, T3_TM.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[a1,a2] = size(A);
[r1,r2,n] = size(X);
[b1,b2] = size(B);
assert(a2 == r1);
assert(r2 == b1);

%% Y(i) = A*X(i)*B
Y = zeros(a1,b2,n);
for i = 1:n
    Y(:,:,i) = A*X(:,:,i)*B;
end

end
