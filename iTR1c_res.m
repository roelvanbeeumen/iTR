function [res,res1,res2,res3,res4] = iTR1c_res(varargin)
%ITR1C_RES  Residual of canonical infinite tensor ring with 1 core.
%
%   res = iTR1c_res(theta,H,Q,S) returns the residual corresponding to the
%   eigenvalue theta and the canonical iTR eigenvector with core Q and
%   singular values matrix S,
%
%       res = res1 + res2 + res3 + res4 - 4*theta*QC
%
%   where
%
%                             [-------]
%                        -----[  Q_C  ]-----
%                      /      [-------]      \
%                      |          |          |
%                    /---\        |          |
%       res1     =  ( H_L )       |          |
%           abc      \---/        |          |
%                      |          |          |
%                      \          c          /
%                        --- a         b ---
%
%                          [-------]       [-------]
%                     -----|  Q_L  |-------|  Q_C  |-----
%                   /      [-------]       [-------]      \
%                   |          |               |          |
%                   |        /-------------------\        |
%       res2     =  |       (          H          )       |
%           abc     |        \-------------------/        |
%                   |          |               |          |
%                   \      [-------]           c          /
%                     -----|  Q_L  |----- a         b ---
%                          [-------]
%
%                          [-------]       [-------]
%                     -----[  Q_C  ]-------[  Q_R  ]-----
%                   /      [-------]       [-------]      \
%                   |          |               |          |
%                   |        /-------------------\        |
%       res3     =  |       (          H          )       |
%           abc     |        \-------------------/        |
%                   |          |               |          |
%                   \          c           [-------]      /
%                     --- a         b -----[  Q_R  ]-----
%                                          [-------]
%
%                          [-------]
%                     -----[  Q_C  ]-----
%                   /      [-------]      \
%                   |          |          |
%                   |          |        /---\
%       res4     =  |          |       ( H_R )
%           abc     |          |        \---/
%                   |          |          |
%                   \          c          /
%                     --- a         b ---
%
%   and the left and right canonical cores are QL = t3_mt(S,Q) and
%   QR = t3_tm(Q,S), respectively, satisfying the following relatations
%
%       t3_tm(QL,S) = t3_mt(S,QR) = t3_mtm(S,Q,S) = QC.
%
%   [res,res1,res2,res3,res4] = iTR1c_res(theta,H,Q,S) also returns the
%   separate terms contributing to the residual.
%
%   res = iTR1c_res(theta,H,QL,C,QR,QC) uses the mixed canonical form where
%   QL, QR, and QC are the left, right, and center canonical cores,
%   respectively. The matrix C is the bound matrix.
%
%   See also iTR1c, iTR1c_HL, iTR1c_HR, iTR1, iTR1_res.

%   Roel Van Beeumen
%   March 18, 2024

if nargin == 4
    % arguments
    theta = varargin{1};
    H = varargin{2};
    Q = varargin{3};
    S = varargin{4};
    % canonical cores
    QL = t3_mt(S,Q);
    QR = t3_tm(Q,S);
    QC = t3_mt(S,QR);
    % residual
    [res,res1,res2,res3,res4] = res_QC(theta,H,QL,S,QR,QC);
elseif nargin == 6
    [res,res1,res2,res3,res4] = res_QC(varargin{:});
else
    error('RVB: only 4 or 5 input arguments are allowed!');
end

end


function [res,res1,res2,res3,res4] = res_QC(theta,H,QL,C,QR,QC)

%% dimensions
[l1,l2,n1] = size(QL);  assert(l1 == l2);
[c1,c2] = size(C);      assert(c1 == c2); assert(c1 == l1);
[r1,r2,n2] = size(QR);  assert(r1 == r2); assert(r1 == l1); assert(n1 == n2);
[q1,q2,n3] = size(QC);  assert(q1 == q2); assert(q1 == l1); assert(n1 == n3);
r = l1;
n = n1;

%% apply H
Z = t4_tt(QL,QC);
HZ = reshape(permute(t4_oper(H,Z),[1,3,2,4]),r*n,r*n);

%% 1st term
HL = iTR1c_HL(H,QL,C);
res1 = t3_mt(HL',QC);

%% 2nd term
HZQL = reshape(permute(QL,[2,1,3]),r,r*n)*HZ;
res2 = reshape(HZQL,r,r,n);

%% 3rd term
HZQR = HZ*reshape(permute(QR,[2,3,1]),r*n,r);
res3 = permute(reshape(HZQR,r,n,r),[1,3,2]);

%% 4th term
HR = iTR1c_HR(H,QR,C);
res4 = t3_tm(QC,HR);

%% residual
res = res1 + res2 + res3 + res4 - 4*theta*QC;

end
