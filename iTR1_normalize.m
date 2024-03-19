function [Xn,Mn] = iTR1_normalize(X,M)
%ITR1_NORMALIZE   Normalize infinite tensor ring with 1 core.
%
%   Xn = iTR1_normalize(X) normalizes the iTR with core X. The core of the
%   normalized iTR is returned in Xn.
%
%   [Xn,Mn] = iTR1_normalize(X,M) normalizes the iTR with core X and matrix
%   M. The core of the normalized iTR is returned in Xn and its matrix in
%   Mn. Note that the Frobenius norm of Mn is 1.
%
%   See also iTR1, iTR1c, iTR2, iTR2c, iTR2_normalize.

%   Roel Van Beeumen
%   March 18, 2024

%% compact
compact = (nargin == 1);

%% dimensions
[x1,x2,~] = size(X);  assert(x1 == x2);
if ~compact
    [m1,m2] = size(M);  assert(m1 == m2); assert(m1 == x1);
end
r = x1;

%% normalize
if compact
    eta = iTR1_tfmat_domeig(X);
    Xn = X/sqrt(eta);
    if nargout > 1
        Mn = eye(r);
    end
else
    eta = iTR1_tfmat_domeig(t3_tm(X,M));
    nrmM = norm(M,'fro');
    Xn = X*(nrmM/sqrt(eta));
    Mn = M/nrmM;
end

end
