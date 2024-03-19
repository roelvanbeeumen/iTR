function [HL,hL] = iTR1c_HL(H,QL,C)
%ITR1C_HL   Inifinte geometric sum of the left transfer matrices.
%
%   HL = iTR1c_HL(H,QL,C) returns the geometric sum of the left transfer
%   matrices
%
%       HL = hL * ( I + TL + TL^2 + ... )
%
%   where
%
%                   [-------]       [-------]
%                ---|  Q_L  |-------|  Q_L  |---
%              /    [-------]       [-------]
%              |        |               |
%              |      /-------------------\
%       hL  =  |     (          H          )
%              |      \-------------------/
%              |        |               |
%              \    [-------]       [-------]
%                ---|  Q_L  |-------|  Q_L  |---
%                   [-------]       [-------]
%
%   and C*C' the dominant right eigenvector of the left transfer matrix.
%
%   [HL,hL] = iTR1c_HL(H,QL,C) also returns the matrix hL.
%
%   See also iTR1c, iTR1c_HR.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[q1,q2,~] = size(QL);  assert(q1 == q2);
[c1,c2] = size(C);     assert(c1 == c2); assert(c1 == q1);
r = q1;

%% left tranfer matrix and its dominant eignevectors
TL = iTR1_tfmat(QL);
vecI = reshape(eye(r),[],1);
vecR = reshape(C*C',[],1);

%% hL
ZL = t4_tt(QL,QL);
HZL = t4_oper(H,ZL);
hL = reshape(permute(HZL,[2,1,3,4]),r,[])*reshape(permute(ZL,[1,3,4,2]),[],r);

%% geometric sum
TLt = TL - vecR*vecI';
HL = reshape(hL(:)'/(eye(r^2) - TLt),r,r);

end
