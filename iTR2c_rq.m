function [theta,theta1,theta2] = iTR2c_rq(H,X,Y,Sxy,Syx)
%ITR2C_RQ   Rayleigh quotient of canonical infinite tensor ring with 2 cores.
%
%   theta = iTR2c_rq(H,X,Y,Sxy,Syx) returns the Rayleigh quotient of the
%   Hamiltonian H and the canonical infinite tensor ring with cores X and Y
%   and singular values matrices Sxy and Syx,
%
%       theta = 1/2*(theta1 + theta2)
%
%   where
%
%                               [-------]           [-------]
%                    --( Syx )--|   X   |--( Sxy )--|   Y   |--( Syx )--
%                  /            [-------]           [-------]            \
%                  |                |                   |                |
%                  |              /-----------------------\              |
%       theta1  =  |             (            H            )             |
%                  |              \-----------------------/              |
%                  |                |                   |                |
%                  \            [-------]           [-------]            /
%                    --( Syx )--|   X   |--( Sxy )--|   Y   |--( Syx )--
%                               [-------]           [-------]
%
%                               [-------]           [-------]
%                    --( Sxy )--|   Y   |--( Syx )--|   X   |--( Sxy )--
%                  /            [-------]           [-------]            \
%                  |                |                   |                |
%                  |              /-----------------------\              |
%       theta2  =  |             (            H            )             |
%                  |              \-----------------------/              |
%                  |                |                   |                |
%                  \            [-------]           [-------]            /
%                    --( Sxy )--|   Y   |--( Syx )--|   X   |--( Sxy )--
%                               [-------]           [-------]
%
%   [theta,theta1,theta2] = iTR2c_rq(H,X,Y,Sxy,Syx) also returns the
%   separate terms contributing to the Rayleigh quotient.
%
%   See also iTR2c, iTR2_rq, iTR1c, iTR1c_rq.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);  assert(x1 == x2);
[y1,y2,n2] = size(Y);  assert(y1 == y2); assert(y1 == x1); assert(n1 == n2);
[s1,s2] = size(Sxy);   assert(s1 == s2); assert(s1 == x1);
[t1,t2] = size(Syx);   assert(t1 == t2); assert(t1 == x1);

%% rayleigh quotient
theta1 = rq(H,X,Y,Sxy,Syx);
theta2 = rq(H,Y,X,Syx,Sxy);
theta = (theta1 + theta2)/2;

end


function theta = rq(H,X,Y,Sxy,Syx)

%% convert to matrix
Z = t4_mtmtm(Syx,X,Sxy,Y,Syx);
Z = reshape(Z,size(Z,1)*size(Z,2),size(Z,3)*size(Z,4));

%% rayleigh quotient
ZtZ = Z'*Z;
theta = sum(sum(ZtZ.*H));

end
