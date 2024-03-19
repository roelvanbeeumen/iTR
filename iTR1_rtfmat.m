function Tr = iTR1_rtfmat(X,S,Y,R)
%ITR1_RTFMAT  Right transfer matrix.
%   ITR1_RTFMAT(X,S) is the right transfer matrix of the infinite tensor ring
%   X-S.
%
%           n                                       x1 [-----] x2   s1       s2
%        -------                                   ----|  X  |---------( S )----
%         \                                            [-----]
%   Tr =   |     [X(:,:,i)*S] (x) [X(:,:,i)*S]  =         | n
%         /                                         x1 [-----] x2   s1       s2
%        -------                                   ----|  X  |---------( S )----
%         i = 1                                        [-----]
%
%   ITR1_RTFMAT(X,S,Y,R) is the right transfer matrix of the infinite tensor
%   rings X-S and Y-R.
%
%           n                                       x1 [-----] x2   s1       s2
%        -------                                   ----|  X  |---------( S )----
%         \                                            [-----]
%   Tr =   |     [Y(:,:,i)*R] (x) [X(:,:,i)*S]  =         | n
%         /                                         y1 [-----] y2   r1       r2
%        -------                                   ----|  Y  |---------( R )----
%         i = 1                                        [-----]
%
%   See also ITR1_TFMAT, ITR1_LTFMAT.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 3, Y = X; end
if nargin < 4, R = S; end

%% dimensions of tensor rings
[x1,x2,n1] = size(X);
[y1,y2,n2] = size(Y);  assert(n1 == n2);
[s1,s2] = size(S);     assert(x2 == s1);
[r1,r2] = size(R);     assert(y2 == r1);

%% transfer matrix
Tr = zeros(x1*y1,s2*r2);
for i = 1:n1
    Tr = Tr + kron(Y(:,:,i)*R,X(:,:,i)*S);
end

end
