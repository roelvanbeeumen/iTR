function [res,resX,resY,resX1,resX2,resX3,resX4,resY1,resY2,resY3,resY4] = ...
    iTR2c_res(theta,H,X,Y,S,T)
%ITR2C_RES  Residual of canonical infinite tensor ring with 2 cores.
%
%   res = iTR2c_res(theta,H,X,Y,S,T) returns the residual corresponding to
%   the eigenvalue theta and the canonical iTR eigenvector with cores X and
%   Y, and singular values matrices S and T,
%
%       res = 1/2*(resX + resY)
%
%   where
%
%       resX = resX1 + resX2 + resX3 + resX4 - 6*theta*TXS
%       resY = resY1 + resY2 + resY3 + resY4 - 6*theta*SYT
%
%   with TXS = t3_mtm(T,X,S), SYT = t3_mtm(S,Y,T)
%
%                                  [-------]
%                         --( T )--[   X   ]--( S )--
%                       /          [-------]          \
%                       |              |              |
%                     /---\            |              |
%       resX1     =  ( HLx )           |              |
%            abc      \---/            |              |
%                       |              |              |
%                       \              c              /
%                         ------- a         b -------
%
%                               [-------]         [-------]
%                      --( S )--|   Y   |--( T )--|   X   |--( S )--
%                    /          [-------]         [-------]          \
%                    |              |                 |              |
%                    |            /---------------------\            |
%       resX2     =  |           (           H           )           |
%            abc     |            \---------------------/            |
%                    |              |                 |              |
%                    \          [-------]             c              /
%                      --( S )--|   Y   |------- a         b -------
%                               [-------]
%
%                               [-------]         [-------]
%                      --( T )--[   X   ]--( S )--[   Y   ]--( T )--
%                    /          [-------]         [-------]          \
%                    |              |                 |              |
%                    |            /---------------------\            |
%       resX3     =  |           (           H           )           |
%            abc     |            \---------------------/            |
%                    |              |                 |              |
%                    \              c             [-------]          /
%                      ------- a         b -------[   Y   ]--( T )--
%                                                 [-------]
%
%                               [-------]
%                      --( T )--[   X   ]--( S )--
%                    /          [-------]          \
%                    |              |              |
%                    |              |            /---\
%       resX4     =  |              |           ( HRx )
%            abc     |              |            \---/
%                    |              |              |
%                    \              c              /
%                      ------- a         b -------
%
%   and
%
%                                  [-------]
%                         --( S )--[   Y   ]--( T )--
%                       /          [-------]          \
%                       |              |              |
%                     /---\            |              |
%       resY1     =  ( HLy )           |              |
%            abc      \---/            |              |
%                       |              |              |
%                       \              c              /
%                         ------- a         b -------
%
%                               [-------]         [-------]
%                      --( T )--|   X   |--( S )--|   Y   |--( T )--
%                    /          [-------]         [-------]          \
%                    |              |                 |              |
%                    |            /---------------------\            |
%       resY2     =  |           (           H           )           |
%            abc     |            \---------------------/            |
%                    |              |                 |              |
%                    \          [-------]             c              /
%                      --( T )--|   X   |------- a         b -------
%                               [-------]
%
%                               [-------]         [-------]
%                      --( S )--[   Y   ]--( T )--[   X   ]--( S )--
%                    /          [-------]         [-------]          \
%                    |              |                 |              |
%                    |            /---------------------\            |
%       resY3     =  |           (           H           )           |
%            abc     |            \---------------------/            |
%                    |              |                 |              |
%                    \              c             [-------]          /
%                      ------- a         b -------[   X   ]--( S )--
%                                                 [-------]
%
%                               [-------]
%                      --( S )--[   Y   ]--( T )--
%                    /          [-------]          \
%                    |              |              |
%                    |              |            /---\
%       resY4     =  |              |           ( HRy )
%            abc     |              |            \---/
%                    |              |              |
%                    \              c              /
%                      ------- a         b -------
%
%   [res,resX,resY,resX1,resX2,resX3,resX4,resY1,resY2,resY3,resY4] = ...
%   iTR2c_res(theta,H,Q,S) also returns the separate terms contributing to
%   the residual.
%
%   See also iTR1c, iTR1c_res, iTR2c_HL, iTR2c_HR.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);  assert(x1 == x2);
[y1,y2,n2] = size(Y);  assert(y1 == y2); assert(y1 == x1); assert(n1 == n2);
[s1,s2] = size(S);     assert(s1 == s2); assert(s1 == x1);
[t1,t2] = size(T);     assert(t1 == t2); assert(t1 == x1);
r = x1;
n = n1;

%% cores
XL = t3_mt(T,X);
XR = t3_tm(X,S);
YL = t3_mt(S,Y);
YR = t3_tm(Y,T);
XC = t3_mtm(T,X,S);
YC = t3_mtm(S,Y,T);

%% apply H
XY = t4_tt(XC,YR);
YX = t4_tt(YL,XC);
HXY = reshape(permute(t4_oper(H,XY),[1,3,2,4]),r*n,r*n);
HYX = reshape(permute(t4_oper(H,YX),[1,3,2,4]),r*n,r*n);

%% 1st terms
[HLx,HLy] = iTR2c_HL(H,X,Y,S,T);
resX1 = t3_mt(HLx',XC);
resY1 = t3_mt(HLy',YC);

%% 2nd terms
YLHYX = reshape(permute(YL,[2,1,3]),r,r*n)*HYX;
resX2 = reshape(YLHYX,r,r,n);
XLHXY = reshape(permute(XL,[2,1,3]),r,r*n)*HXY;
resY2 = reshape(XLHXY,r,r,n);

%% 3rd terms
YRHXY = HXY*reshape(permute(YR,[2,3,1]),r*n,r);
resX3 = permute(reshape(YRHXY,r,n,r),[1,3,2]);
XRHYX = HYX*reshape(permute(XR,[2,3,1]),r*n,r);
resY3 = permute(reshape(XRHYX,r,n,r),[1,3,2]);

%% 4th terms
[HRx,HRy] = iTR2c_HR(H,X,Y,S,T);
resX4 = t3_tm(XC,HRx);
resY4 = t3_tm(YC,HRy);

%% residual
resX = resX1 + resX2 + resX3 + resX4 - 6*theta*XC;
resY = resY1 + resY2 + resY3 + resY4 - 6*theta*YC;
res = (resX + resY)/2;

end
