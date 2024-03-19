function Tl = iTR1_ltfmat(S,X,R,Y)
%ITR1_LTFMAT  Left transfer matrix.
%   ITR1_LTFMAT(S,X) is the left transfer matrix of the infinite tensor ring
%   S-X.
%
%           n                                       s1       s2   x1 [-----] x2
%        -------                                   ----( S )---------|  X  |----
%         \                                                          [-----]
%   Tl =   |     [S*X(:,:,i)] (x) [S*X(:,:,i)]  =                       | n
%         /                                         s1       s2   x1 [-----] x2
%        -------                                   ----( S )---------|  X  |----
%         i = 1                                                      [-----]
%
%   ITR1_LTFMAT(S,X,R,Y) is the left transfer matrix of the infinite tensor
%   rings S-X and R-Y.
%
%           n                                       s1       s2   x1 [-----] x2
%        -------                                   ----( S )---------|  X  |----
%         \                                                          [-----]
%   Tl =   |     [R*Y(:,:,i)] (x) [S*X(:,:,i)]  =                       | n
%         /                                         r1       r2   y1 [-----] y2
%        -------                                   ----( R )---------|  Y  |----
%         i = 1                                                      [-----]
%
%   See also ITR1_TFMAT, ITR1_RTFMAT.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 3, R = S; end
if nargin < 4, Y = X; end

%% dimensions of tensor rings
[x1,x2,n1] = size(X);
[y1,y2,n2] = size(Y);  assert(n1 == n2);
[s1,s2] = size(S);     assert(s2 == x1);
[r1,r2] = size(R);     assert(r2 == y1);

%% transfer matrix
Tl = zeros(s1*r1,x2*y2);
for i = 1:n1
    Tl = Tl + kron(R*Y(:,:,i),S*X(:,:,i));
end

end
