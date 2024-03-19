function [Xn,Yn,Cn,Dn] = iTR2_normalize(X,Y,C,D)
%ITR2_NORMALIZE   Normalize infinite tensor ring with 2 cores.
%
%   [Xn,Yn] = iTR2_normalize(X,Y) normalizes the iTR with cores X and Y.
%   The cores of the normalized iTR are returned in Xn and Yn,
%   respectively.
%
%   [Xn,Yn,Cn,Dn] = iTR2_normalize(X,Y,C,D) normalizes the iTR with cores X
%   and Y, and matrices C and D. The cores of the normalized iTR are
%   returned in Xn and Yn, respectively, and its matrices in Cn and Dn,
%   respectively. Note that the Frobenius norms of Cn and Dn are 1.
%
%   See also iTR2, iTR2c, iTR1, iTR1c, iTR1_normalize.

%   Roel Van Beeumen
%   March 18, 2024

%% compact
compact = (nargin == 2);

%% dimensions
[x1,x2,n1] = size(X);   assert(x1 == x2);
[y1,y2,n2] = size(Y);   assert(y1 == y2); assert(y1 == x1); assert(n1 == n2);
if ~compact
    [c1,c2] = size(C);  assert(c1 == c2); assert(c1 == x1);
    [d1,d2] = size(D);  assert(d1 == d2); assert(d1 == x1);
end
r = x1;

%% iTR1
if compact
    Z = t3_tt(X,Y);
else
    Z = t3_tt(t3_tm(X,C),t3_tm(Y,D));
end
eta = iTR1_tfmat_domeig(Z);

%% normalize
if compact
    scal = sqrt(sqrt(eta));
    Xn = X/scal;
    Yn = Y/scal;
    if nargout > 2
        Cn = eye(r);
    elseif nargout > 3
        Dn = eye(r);
    end
else
    nrmC = norm(C,'fro');
    nrmD = norm(D,'fro');
    scal = sqrt(nrmC)*sqrt(nrmD)/sqrt(sqrt(eta));
    Xn = X*scal;
    Yn = Y*scal;
    Cn = C/nrmC;
    Dn = D/nrmD;
end

end
