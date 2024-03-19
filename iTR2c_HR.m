function [HRx,HRy,hRxy,hRyx] = iTR2c_HR(H,X,Y,S,T)
%ITR2C_HR   Inifinte geometric sums of the right transfer matrices.
%
%   [HRx,HRy] = iTR2c_HR(H,X,Y,S,T) returns the geometric sums of the right
%   transfer matrices
%
%       HRx = ( I + TR_YX^2 + ... ) * (hRyx + TR_Y*hRxy)
%       HRy = ( I + TR_XY^2 + ... ) * (hRxy + TR_X*hRyx)
%
%   where
%
%                   [-------]         [-------]
%                ---|   X   |--( S )--|   Y   |--( T )--
%                   [-------]         [-------]          \
%                       |                 |              |
%                     /---------------------\            |
%       hRxy  =      (           H           )           |
%                     \---------------------/            |
%                       |                 |              |
%                   [-------]         [-------]          /
%                ---|   X   |--( S )--|   Y   |--( T )--
%                   [-------]         [-------]
%
%                   [-------]         [-------]
%                ---|   Y   |--( T )--|   X   |--( S )--
%                   [-------]         [-------]          \
%                       |                 |              |
%                     /---------------------\            |
%       hRyx  =      (           H           )           |
%                     \---------------------/            |
%                       |                 |              |
%                   [-------]         [-------]          /
%                ---|   Y   |--( T )--|   X   |--( S )--
%                   [-------]         [-------]
%
%   with T*T' and S'*S the dominant left eigenvectors of the right transfer
%   matrices TR_XY and TR_YX, respectively.
%
%   [HRx,HRy,hRxy,hRyx] = iTR2c_HR(H,X,Y,S,T) also returns the matrices
%   hRxy and hRyx.
%
%   See also iTR2c, iTR2c_HL, iTR1c, iTR1c_HL, iTR1c_HR.

%   Roel Van Beeumen
%   March 18, 2024

%% dimensions
[x1,x2,n1] = size(X);  assert(x1 == x2);
[y1,y2,n2] = size(Y);  assert(y1 == y2); assert(y1 == x1); assert(n1 == n2);
[s1,s2] = size(S);     assert(s1 == s2); assert(s1 == x1);
[t1,t2] = size(T);     assert(t1 == t2); assert(t1 == x1);
r = x1;

%% right cores
XR = t3_tm(X,S);
YR = t3_tm(Y,T);

%% hRxy and hRyx
XY = t4_tt(XR,YR);
YX = t4_tt(YR,XR);
HXY = t4_oper(H,XY);
HYX = t4_oper(H,YX);
hRxy = reshape(HXY,r,[])*reshape(permute(XY,[2,3,4,1]),[],r);
hRyx = reshape(HYX,r,[])*reshape(permute(YX,[2,3,4,1]),[],r);

%% vectorize
vecI = reshape(eye(r),[],1);
vecS2 = reshape(S^2,[],1);
vecT2 = reshape(T^2,[],1);
hRxyt = hRxy(:);
hRyxt = hRyx(:);

%% geometric sums
% HRx = reshape((eye(r^2) - TYXt)\(hRyx + TY*hRxy),r,r);
HRx = reshape(mymldivide(...
    @(x) x - iTR1_tfmatvec(YR,iTR1_tfmatvec(XR,x)) + vecI*(vecS2'*x),...
    hRyxt + iTR1_tfmatvec(YR,hRxyt)),r,r);
% HRy = reshape((eye(r^2) - TXYt)\(hRxy + TX*hRyx),r,r);
HRy = reshape(mymldivide(...
    @(x) x - iTR1_tfmatvec(XR,iTR1_tfmatvec(YR,x)) + vecI*(vecT2'*x),...
    hRxyt + iTR1_tfmatvec(XR,hRyxt)),r,r);

end


function C = mymldivide(Afun,B)

[C,flag] = gmres(Afun,B,10,1e-8,100);
if flag
    warning('RVB: gmres did not converge!');
end

end


function y = iTR1_tfmatvec(X,v)

%% dimensions
[r,~,n] = size(X);

%% reshape
M = reshape(v,r,r);

%% matvec
y = zeros(size(v));
for i = 1:n
    Y = X(:,:,i)*M*X(:,:,i)';
    y = y + Y(:);
end

end
