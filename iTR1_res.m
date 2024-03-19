function [res,res1,res2,res3,res4] = iTR1_res(theta,H,X,C,VR,WR,VL,WL)
%ITR1_RES   Residual of infinite tensor ring with 1 core.
%
%   res = iTR1_res(theta,H,X,C) returns the residual corresponding to the
%   eigenvalue theta and the iTR eigenvector with core X and matrix C,
%
%       res = 1/2*(res1 + res2 + res3 + res4) - theta*CXC
%
%   where CXC = t3_mtm(C,X,C),
%
%                                 [-------]         [-------]
%                        --( C )--|   X   |--( C )--|   X   |--( C )--
%                      /          [-------]         [-------]          \
%                      |              |                 |              |
%                    /---\          /---------------------\          /---\
%       res1     =  ( W_L )        (           H           )        ( V_R )
%           abc      \---/          \---------------------/          \---/
%                      |              |                 |              |
%                      \          [-------]             c              /
%                        --( C )--|   X   |------- a         b -------
%                                 [-------]
%
%                                 [-------]         [-------]
%                        --( C )--[   X   ]--( C )--[   X   ]--( C )--
%                      /          [-------]         [-------]          \
%                      |              |                 |              |
%                    /---\          /---------------------\          /---\
%       res2     =  ( W_L )        (           H           )        ( V_R )
%           abc      \---/          \---------------------/          \---/
%                      |              |                 |              |
%                      \              c             [-------]          /
%                        ------- a         b -------[   X   ]--( C )--
%                                                   [-------]
%
%                                 [-------]
%                        --( C )--[   X   ]--( C )--
%                      /          [-------]          \
%                      |              |              |
%                    /---\            |            /---\
%       res3     =  ( H_L )           |           ( V_R )
%           abc      \---/            |            \---/
%                      |              |              |
%                      \              c              /
%                        ------- a         b -------
%
%                                 [-------]
%                        --( C )--[   X   ]--( C )--
%                      /          [-------]          \
%                      |              |              |
%                    /---\            |            /---\
%       res4     =  ( W_L )           |           ( H_R )
%           abc      \---/            |            \---/
%                      |              |              |
%                      \              c              /
%                        ------- a         b -------
%
%   WL the left eigenvector corresponding to the dominant eigenvalue 1 of
%   the left transfer matrix, and VR the right eigenvector corresponing to
%   the dominant eigenvalue 1 of the right transfer matrix.
%
%   [res,res1,res2,res3,res4] = iTR1_res(theta,H,X,C) also returns the
%   separate terms contributing to the residual.
%
%   iTR1_res(theta,H,X,C,VR,WR,VL,WL) uses the given matrices WR and VR as
%   left and right eigenvectors, respectively, corresponding to the
%   dominant eigenvalue 1 of the right transfer matrix, and the given
%   matrices WL and VL as left and right eigenvectors, repsectively,
%   corresponding to the dominant eigenvalue 1 of the left transfer matrix.
%   Note that VR and WL should be normalized so that trace(WR'*VR) = 1,
%   trace(WL'*VL) = 1, and trace(WL'*C*VR*C') = 1.
%
%   See also iTR1, iTR1c_HL, iTR1c_HR, iTR1c, iTR1c_res.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n] = size(X);       assert(x1 == x2);
[c1,c2] = size(C);         assert(c1 == c2); assert(c1 == x1);
if nargin == 8
    [vr1,vr2] = size(VR);  assert(vr1 == vr2); assert(vr1 == x1);
    [wr1,wr2] = size(WR);  assert(wr1 == wr2); assert(wr1 == x1);
    [vl1,vl2] = size(VL);  assert(vl1 == vl2); assert(vl1 == x1);
    [wl1,wl2] = size(WL);  assert(wl1 == wl2); assert(wl1 == x1);
end
r = x1;

%% left, right, and center cores
XL = t3_mt(C,X);
XR = t3_tm(X,C);
XC = t3_mt(C,XR);

%% default arguments
if nargin < 8
    % dominant eigenvectors of left transfer matrix
    [eta,VL,WL] = iTR1_tfmat_domeig(XL,'leftnormalize');
    if abs(eta - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
    % dominant eigenvectors of right transfer matrix
    [eta,VR,WR] = iTR1_tfmat_domeig(XR,'rightnormalize');
    if abs(eta - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
    % normalize so that trace(WL'*C*VR*C') = 1
    scal = trace(WL'*C*VR*C');
    WL = WL/scal;
    VL = VL*scal;
end

%% apply H
Z = t4_mtmtm(WL',XL,eye(r),XC,VR);
HZ = reshape(permute(t4_oper(H,Z),[1,3,2,4]),r*n,r*n);

%% 1st term
HZXL = reshape(permute(XL,[2,1,3]),r,r*n)*HZ;
res1 = reshape(HZXL,r,r,n);

%% 2nd term
HZXR = HZ*reshape(permute(XR,[2,3,1]),r*n,r);
res2 = permute(reshape(HZXR,r,n,r),[1,3,2]);

%% 3rd term
HL = iTR1_HL(H,X,C,VL,WL);
res3 = t3_mtm(HL',XC,VR);

%% 4th term
HR = iTR1_HR(H,X,C,VR,WR);
res4 = t3_mtm(WL',XC,HR);

%% residual
res = (res1 + res2 + res3 + res4)/2 - theta*XC;

end
