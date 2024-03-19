function [nrm,nrm1,nrm2] = iTR2c_norm(X,Y,Sxy,Syx)
%ITR2C_NORM   Norm of canonical infinite tensor ring with 2 cores.
%
%   nrm = iTR2c_norm(X,Y,Sxy,Syx) returns the norm of the canonical iTR
%   with cores X and Y, and singular values matrices Sxy and Syx,
%
%       nrm = 1/2*(nrm1 + nrm2)
%
%   where
%
%                             [-------]           [-------]
%                  --( Syx )--|   X   |--( Sxy )--|   Y   |--( Syx )--
%                /            [-------]           [-------]            \
%       nrm1  =  |                |                   |                |
%                \            [-------]           [-------]            /
%                  --( Syx )--|   X   |--( Sxy )--|   Y   |--( Syx )--
%                             [-------]           [-------]
%
%                             [-------]           [-------]
%                  --( Sxy )--|   Y   |--( Syx )--|   X   |--( Sxy )--
%                /            [-------]           [-------]            \
%       nrm2  =  |                |                   |                |
%                \            [-------]           [-------]            /
%                  --( Sxy )--|   Y   |--( Syx )--|   X   |--( Sxy )--
%                             [-------]           [-------]
%
%   [nrm,nrm1,nrm2] = iTR2c_norm(X,Y,Sxy,Syx) also returns the separate
%   terms contributing to the norm.
%
%   See also iTR2c, iTR2c_rq, iTR1c, iTR1c_norm, iTR1c_rq.

%   Roel Van Beeumen
%   March 18, 2024

%% norm
nrm1 = iTR2norm(X,Y,Sxy,Syx);
nrm2 = iTR2norm(Y,X,Syx,Sxy);
nrm = (nrm1 + nrm2)/2;

end


function nrm = iTR2norm(X,Y,Sxy,Syx)

%% convert to matrix
Z = t4_mtmtm(Syx,X,Sxy,Y,Syx);
Z = reshape(Z,size(Z,1)*size(Z,2),size(Z,3)*size(Z,4));

%% norm
nrm = trace(Z'*Z);

end
