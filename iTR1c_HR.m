function [HR,hR] = iTR1c_HR(H,QR,C)
%ITR1C_HR   Inifinte geometric sum of the right transfer matrices.
%
%   HR = iTR1c_HR(H,QR,C) returns the geometric sum of the right transfer
%   matrices
%
%       HR = ( I + TR + TR^2 + ... ) * hR
%
%   where
%
%                 [-------]       [-------]
%              ---|  Q_R  |-------|  Q_R  |---
%                 [-------]       [-------]    \
%                     |               |        |
%                   /-------------------\      |
%       hR  =      (          H          )     |
%                   \-------------------/      |
%                     |               |        |
%                 [-------]       [-------]    /
%              ---|  Q_R  |-------|  Q_R  |---
%                 [-------]       [-------]
%
%   and C'*C the dominant left eigenvector of the right transfer matrix.
%
%   [HR,hR] = iTR1c_HR(H,QR,C) also returns the matrix hR.
%
%   See also iTR1c, iTR1c_HL.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[q1,q2,~] = size(QR);  assert(q1 == q2);
[c1,c2] = size(C);     assert(c1 == c2); assert(c1 == q1);
r = q1;

%% right tranfer matrix and its dominant eignevectors
TR = iTR1_tfmat(QR);
vecI = reshape(eye(r),[],1);
vecL = reshape(C'*C,[],1);

%% hR
ZR = t4_tt(QR,QR);
HZR = t4_oper(H,ZR);
hR = reshape(HZR,r,[])*reshape(permute(ZR,[2,3,4,1]),[],r);

%% geometric sum
TRt = TR - vecI*vecL';
HR = reshape((eye(r^2) - TRt)\hR(:),r,r);

end
