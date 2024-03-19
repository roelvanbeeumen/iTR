function [theta,theta1,theta2] = iTR2_rq(H,X,Y,C,D,Vxy,Wxy,Vyx,Wyx)
%ITR2_RQ  Rayleigh quotient of infinite tensor ring with 2 cores.
%
%   theta = iTR2_rq(H,X,Y,C,D) returns the Rayleigh quotient of the
%   Hamiltonian H and the canonical infinite tensor ring with cores X and Y
%   and matrices C and D,
%
%       theta = 1/2*(theta1 + theta2)
%
%   where
%
%                                [-------]         [-------]
%                       --( D )--|   X   |--( C )--|   Y   |--( D )--
%                     /          [-------]         [-------]          \
%                     |              |                 |              |
%                   /---\          /---------------------\          /---\
%       theta1  =  ( Wxy )        (           H           )        ( Vxy )
%                   \---/          \---------------------/          \---/
%                     |              |                 |              |
%                     \          [-------]         [-------]          /
%                       --( D )--|   X   |--( C )--|   Y   |--( D )--
%                                [-------]         [-------]
%
%                                [-------]         [-------]
%                       --( C )--|   Y   |--( D )--|   X   |--( C )--
%                     /          [-------]         [-------]          \
%                     |              |                 |              |
%                   /---\          /---------------------\          /---\
%       theta2  =  ( Wyx )        (           H           )        ( Vyx )
%                   \---/          \---------------------/          \---/
%                     |              |                 |              |
%                     \          [-------]         [-------]          /
%                       --( C )--|   Y   |--( D )--|   X   |--( C )--
%                                [-------]         [-------]
%
%   [theta,theta1,theta2] = iTR2_rq(H,X,Y,C,D) also returns the separate
%   terms contributing to the Rayleigh quotient.
%
%   iTR2_rq(H,X,Y,C,D,Vxy,Wxy,Vyx,Wyx) uses the given matrix Vxy as the
%   right eigenvector corrsponding to the dominant eigenvalue 1 of the
%   right transfer matrix T_XCYD, the matrix Wxy as the left eigenvector
%   corresponding to the dominant eigenvalue 1 of the left transfer matrix
%   T_DXCY, the matrix Vyx as the right eigenvector corresponding to the
%   dominant eigenvalue 1 of the right transfer matrix T_YDXC, and the
%   matrix Wyx as the left eigenvector corresponding to the dominant
%   eigenvalue 1 of the left transfer matrix T_CYDX. Note that Vxy and Wxy
%   should be normalized so that trace(Wxy'*D*Vxy*D') = 1, and Vyx and Wyx
%   so that trace(Wyx'*C*Vyx*C') = 1.
%
%   See also iTR2, iTR2c, iTR2c_rq, iTR1, iTR1_rq, iTR1c, iTR1c_rq.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);  assert(x1 == x2);
[y1,y2,n2] = size(Y);  assert(y1 == y2); assert(y1 == x1); assert(n1 == n2);
[c1,c2] = size(C);     assert(c1 == c2); assert(c1 == x1);
[d1,d2] = size(D);     assert(d1 == d2); assert(d1 == x1);

%% rayleigh quotient
if nargin < 9
    theta1 = rq(H,X,Y,C,D);
    theta2 = rq(H,Y,X,D,C);
else
    theta1 = rq(H,X,Y,C,D,Vxy,Wxy);
    theta2 = rq(H,Y,X,D,C,Vyx,Wyx);
end
theta = (theta1 + theta2)/2;

end


function theta = rq(H,X,Y,C,D,VR,WL)

%% dimension
r = size(X,1);

%% left and right cores
ZL = t3_mtmt(D,X,C,Y);  % D-X-C-Y
ZR = t3_tmtm(X,C,Y,D);  % X-C-Y-D

%% default arguments
if nargin < 7
    % left dominant eigenvector
    [etaL,WL] = iTR1_tfmat_domeig(ZL,'L');
    if abs(etaL - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
    % right dominant eigenvector
    [etaR,VR] = iTR1_tfmat_domeig(ZR,'R');
    if abs(etaR - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
    % normalize so that trace(WL'*D*VR*D') = 1
    WL = WL/trace(WL'*D*VR*D');
end

%% convert to matrices
DXCYDVR = reshape(t3_tm(ZL,D*VR),r^2,[]);
WLDXCYD = reshape(t3_mt(WL*D,ZR),r^2,[]);

%% raleigh quotient
theta = sum(sum((WLDXCYD'*DXCYDVR).*H));

end
