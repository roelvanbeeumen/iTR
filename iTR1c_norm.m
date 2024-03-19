function nrm = iTR1c_norm(Q,S)
%ITR1C_NORM   Norm of canonical infinite tensor ring with 1 core.
%
%   nrm = iTR1c_norm(Q,S) returns the norm of the canonical iTR with core Q
%   and singular values matrix S,
%
%                           [-------]
%                 --( S )---|   Q   |---( S )--
%               /           [-------]           \
%       nrm  =  |               |               |
%               \           [-------]           /
%                 --( S )---|   Q   |---( S )--
%                           [-------]
%
%   See also iTR1c, iTR1c_rq, iTR2c, iTR2c_norm, iTR2c_rq.

%   Roel Van Beeumen
%   March 18, 2024

%% convert to matrix
Z = t3_mtm(S,Q,S);
Z = reshape(Z,size(Z,1)*size(Z,2),size(Z,3));

%% norm
nrm = trace(Z'*Z);

end
