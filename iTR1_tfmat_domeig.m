function [eta,V,W] = iTR1_tfmat_domeig(X,type)
%ITR1_TFMAT_DOMEIG  Dominant eigenpair of transfer matrix.
%
%   eta = iTR1_tfmat_domeig(X) returns the dominant eigenvalue of the
%   transfer matrix
%
%              n                               x1 [-----] x2
%           -------                           ----|  X  |----
%            \                                    [-----]
%       T =   |     X(:,:,i) (x) X(:,:,i)  =         | n       .
%            /                                 x1 [-----] x2
%           -------                           ----|  X  |----
%            i = 1                                [-----]
%
%   [eta,V,W] = iTR1_tfmat_domeig(X) returns the dominant eigenvalue eta
%   and corresponding left and right eigenvectors, W and V, respectively,
%   of the transfer matrix T. The eigenvectors are normalized such that
%   trace(W'*V) = 1.
%
%   [eta,V,W] = iTR1_tfmat_domeig(X,NORMALIZATION) normalizes the left and
%   right eigenvectors according to
%
%            'default' - normalization: trace(W'*V) = 1
%      'leftnormalize' - normalization: trace(W'*V) = 1 and W(1,1) = 1
%     'rightnormalize' - normalization: trace(W'*V) = 1 and V(1,1) = 1
%
%   [eta,W] = iTR1_tfmat_domeig(X,'L') only returns the dominant eigenvalue
%   eta and corresponding left eigenvector W.
%
%   [eta,V] = iTR1_tfmat_domeig(X,'R') only returns the dominant eigenvalue
%   eta and corresponding right eigenvector V.
%
%   See also iTR1_tfmat.

%   Roel Van Beeumen
%   March 18, 2024

%% default argument
if nargin < 2, type = 'none'; end

%% domeig
if nargout < 2
    eta = domeig(X);
elseif strcmpi(type,'l')
    [eta,V] = domeig_left(X);
elseif strcmpi(type,'r')
    [eta,V] = domeig_right(X);
else
    [eta,V,W] = domeig_both(X,type);
end

%% check for real eta
assert(isreal(eta) || abs(imag(eta)) < 100*eps);
eta = real(eta);

end


function eta = domeig(X)

%% rank
r = size(X,1);

%% small
small = r^2 < 100;

%% dominant eigenpair
if small
    T = iTR1_tfmat(X);
    d = eig(T);
    [~,i] = max(abs(d));
    eta = d(i);
else
    eta = eigs(@(v) matvec(X,v),r^2,1,'largestreal');
end

end


function [eta,V,W] = domeig_both(X,normalization)

%% rank
r = size(X,1);

%% small
small = r^2 < 100;

%% dominant eigenpair
if small
    T = iTR1_tfmat(X);
    [V,d,W] = eig(T,'vector');
    [~,i] = max(abs(d));
    V = reshape(V(:,i),r,r);
    W = reshape(W(:,i),r,r);
    eta = d(i);
else
    [V,dr] = eigs(@(v) matvec(X,v),r^2,1,'largestreal');
    [W,dl] = eigs(@(v) matvect(X,v),r^2,1,'largestreal');
    V = reshape(V,r,r);
    W = reshape(W,r,r);
    eta = (dr + dl)/2;
end
V = (V + V')/2;
W = (W + W')/2;

%% normalization
if strcmpi(normalization,'leftnormalize')
    W = W/W(1,1);
    V = V/trace(W'*V);
elseif strcmpi(normalization,'rightnormalize')
    V = V/V(1,1);
    W = W/trace(W'*V);
else
    W = W/trace(W'*V);
end

end


function [eta,W] = domeig_left(X)

%% rank
r = size(X,1);

%% small
small = r^2 < 100;

%% dominant eigenpair
if small
    T = iTR1_tfmat(X);
    [W,d] = eig(T','vector');
    [~,i] = max(abs(d));
    W = reshape(W(:,i),r,r);
    eta = d(i);
else
    [W,eta] = eigs(@(v) matvect(X,v),r^2,1,'largestreal');
    W = reshape(W,r,r);
end
W = (W + W')/2;

end


function [eta,V] = domeig_right(X)

%% rank
r = size(X,1);

%% small
small = r^2 < 100;

%% dominant eigenpair
if small
    T = iTR1_tfmat(X);
    [V,d] = eig(T,'vector');
    [~,i] = max(abs(d));
    V = reshape(V(:,i),r,r);
    eta = d(i);
else
    [V,eta] = eigs(@(v) matvec(X,v),r^2,1,'largestreal');
    V = reshape(V,r,r);
end
V = (V + V')/2;

end


function y = matvec(X,v)

%% dimensions
[r,~,n] = size(X);

%% reshape
M = reshape(v,r,r);

%% matvec
y = zeros(size(v));
for i = 1:n
    Y = X(:,:,i)*M*X(:,:,i)';
    y = y + Y(:);
end

end


function y = matvect(X,v)

%% dimensions
[r,~,n] = size(X);

%% reshape
M = reshape(v,r,r);

%% matvec
y = zeros(size(v));
for i = 1:n
    Y = X(:,:,i)'*M*X(:,:,i);
    y = y + Y(:);
end

end
