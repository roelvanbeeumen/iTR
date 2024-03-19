function T = iTR1_tfmat(X,Y)
%ITR1_TFMAT   Transfer matrix.
%   ITR1_TFMAT(X) is the transfer matrix of the infinite tensor ring X.
%
%          n                               x1 [-----] x2
%       -------                           ----|  X  |----
%        \                                    [-----]
%   T =   |     X(:,:,i) (x) X(:,:,i)  =         | n
%        /                                 x1 [-----] x2
%       -------                           ----|  X  |----
%        i = 1                                [-----]
%
%   ITR1_TFMAT(X,Y) is the transfer matrix of the infinite tensor rings X and Y.
%
%          n                               x1 [-----] x2
%       -------                           ----|  X  |----
%        \                                    [-----]
%   T =   |     Y(:,:,i) (x) X(:,:,i)  =         | n
%        /                                 y1 [-----] y2
%       -------                           ----|  Y  |----
%        i = 1                                [-----]
%
%   See also ITR1_LTFMAT, ITR1_RTFMAT.

%   Roel Van Beeumen
%   March 18, 2024

%% default argument
if nargin < 2, Y = X; end

%% dimensions of tensor rings
[x1,x2,n1] = size(X);
[y1,y2,n2] = size(Y);
assert(n1 == n2);

%% transfer matrix
T = zeros(x1*y1,x2*y2);
for i = 1:n1
    T = T + kron(Y(:,:,i),X(:,:,i));
end

end
