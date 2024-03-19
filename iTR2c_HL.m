function [HLx,HLy,hLxy,hLyx] = iTR2c_HL(H,X,Y,S,T)
%ITR2C_HL   Inifinte geometric sums of the left transfer matrices.
%
%   [HLx,HLy] = iTR2c_HL(H,X,Y,S,T) returns the geometric sums of the left
%   transfer matrices
%
%       HLx = (hLxy + hLyx*TL_Y) * ( I + TL_XY^2 + ... )
%       HLy = (hLyx + hLxy*TL_X) * ( I + TL_YX^2 + ... )
%
%   where
%
%                           [-------]         [-------]
%                  --( T )--|   X   |--( S )--|   Y   |---
%                /          [-------]         [-------]
%                |              |                 |
%                |            /---------------------\
%       hLxy  =  |           (           H           )
%                |            \---------------------/
%                |              |                 |
%                \          [-------]         [-------]
%                  --( T )--|   X   |--( S )--|   Y   |---
%                           [-------]         [-------]
%
%                           [-------]         [-------]
%                  --( S )--|   Y   |--( T )--|   X   |---
%                /          [-------]         [-------]
%                |              |                 |
%                |            /---------------------\
%       hLyx  =  |           (           H           )
%                |            \---------------------/
%                |              |                 |
%                \          [-------]         [-------]
%                  --( S )--|   Y   |--( T )--|   X   |---
%                           [-------]         [-------]
%
%   with T*T' and S*S' the dominant right eigenvectors of the left transfer
%   matrices TL_XY and TL_YX, respectively.
%
%   [HLx,HLy,hLxy,hLyx] = iTR2c_HL(H,X,Y,S,T) also returns the matrices
%   hLxy and hLyx.
%
%   See also iTR2c, iTR2c_HR, iTR1c, iTR1c_HL, iTR1c_HR.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);  assert(x1 == x2);
[y1,y2,n2] = size(Y);  assert(y1 == y2); assert(y1 == x1); assert(n1 == n2);
[s1,s2] = size(S);     assert(s1 == s2); assert(s1 == x1);
[t1,t2] = size(T);     assert(t1 == t2); assert(t1 == x1);
r = x1;

%% left cores
XL = t3_mt(T,X);
YL = t3_mt(S,Y);

%% hLxy and hLyx
XY = t4_tt(XL,YL);
YX = t4_tt(YL,XL);
HXY = t4_oper(H,XY);
HYX = t4_oper(H,YX);
hLxy = reshape(permute(HXY,[2,1,3,4]),r,[])*reshape(permute(XY,[1,3,4,2]),[],r);
hLyx = reshape(permute(HYX,[2,1,3,4]),r,[])*reshape(permute(YX,[1,3,4,2]),[],r);

%% vectorize
vecI = reshape(eye(r),[],1);
vecS2 = reshape(S^2,[],1);
vecT2 = reshape(T^2,[],1);

%% geometric sums
% HLx = reshape((hLxy + hLyx*TY)/(eye(r^2) - TXYt),r,r);
HLx = reshape(mymrdivide(...
    hLxy(:) + iTR1_tfmatvect(YL,hLyx(:)),...
    @(x) x - iTR1_tfmatvect(YL,iTR1_tfmatvect(XL,x)) + vecI*(vecT2'*x)),r,r);
% HLy = reshape((hLyx + hLxy*TX)/(eye(r^2) - TYXt),r,r);
HLy = reshape(mymrdivide(...
    hLyx(:) + iTR1_tfmatvect(XL,hLxy(:)),...
    @(x) x - iTR1_tfmatvect(XL,iTR1_tfmatvect(YL,x)) + vecI*(vecS2'*x)),r,r);

end


function C = mymrdivide(AT,BTfun)
% AT = A.'
% BTfun = @(x) B.'*x

[C,flag] = gmres(@(x) BTfun(x),AT,10,1e-8,100);
C = C.';
if flag
    warning('RVB: gmres did not converge!');
end

end


function y = iTR1_tfmatvect(X,v)

%% dimensions
[r,~,n] = size(X);

%% reshape
M = reshape(v,r,r);

%% matvec
y = zeros(size(v));
for i = 1:n
    Y = X(:,:,i)'*M*X(:,:,i);
    y = y + Y(:);
end

end
