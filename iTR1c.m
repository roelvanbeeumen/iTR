function [Q,S] = iTR1c(A,M)
%ITR1C  Canonical form of infinite tensor ring with 1 core.
%
%               [-------]
%     --( S )---|   Q   |---( S )--
%   /           [-------]           \
%   |               |               |  =  1
%   \           [-------]           /
%     --( S )---|   Q   |---( S )--
%               [-------]
%
%   See also ITR1, ITR2, ITR2C.

%   Reference:
%   R. Orus and G. Vidal. Infinite time-evolving block decimation algorithm
%   beyond unitary evolution. Physical Review B 78, 155117, 2008.
%
%   Roel Van Beeumen
%   March 18, 2024

%% compact
compact = (nargin == 1);

%% dimensions
[r1,r2,~] = size(A);
if compact
    assert(r1 == r2);
else
    [m1,m2] = size(M);
    assert(r1 == m2);
    assert(r2 == m1);
end

%% dominant eigenvectors of transfer matrix/matrices
if compact
    [eta,Vr,Vl] = iTR1_tfmat_domeig(A);
else
    [etaR,Vr] = iTR1_tfmat_domeig(t3_tm(A,M),'R');
    [etaL,Vl] = iTR1_tfmat_domeig(t3_mt(M,A),'L');
    eta = (etaR + etaL)/2;
end

%% eigenvalue decompositions of Vl and Vr
[Ur,lamr] = eig(Vr,'vector'); lamr = abs(lamr);
[Ul,laml] = eig(Vl,'vector'); laml = abs(laml);

%% singular value decomposition
if compact
    LUUL = (diag(sqrt(laml))*Ul')*(Ur*diag(sqrt(lamr)));
else
    LUUL = (diag(sqrt(laml))*Ul')*M*(Ur*diag(sqrt(lamr)));
end
[V,S,W] = svd(LUUL,'econ');
nrmS = norm(S,'fro');

%% canonical form
Sl = W'*diag(1./sqrt(lamr))*Ur';
Sr = Ul*diag(1./sqrt(laml))*V;
Q = t3_mtm(Sl,A,Sr);

%% normalize
Q = Q*(nrmS/sqrt(eta));
S = S/nrmS;

end
