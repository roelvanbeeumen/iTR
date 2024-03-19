function theta = iTR1_rq(varargin)
%ITR1_RQ  Rayleigh quotient of infinite tensor ring with 1 core.
%
%   iTR1_rq(H,X) returns the Rayleigh quotient
%
%                           [-------]       [-------]
%                      -----|   X   |-------|   X   |-----
%                    /      [-------]       [-------]      \
%                    |          |               |          |
%                  /---\      /-------------------\      /---\
%       theta  =  (  W  )    (          H          )    (  V  )
%                  \---/      \-------------------/      \---/
%                    |          |               |          |
%                    \      [-------]       [-------]      /
%                      -----|   X   |-------[   X   ]-----
%                           [-------]       [-------]
%
%   where W and V are the left and right eigenvectors, respectively,
%   corresponding to the dominant eigenvalue 1 of the transfer matrix.
%
%   iTR1_rq(H,X,V,W) uses the given matrices V and W. Note that V and W
%   should be normalized so that trace(W'*V) = 1.
%
%   iTR1_rq(H,X,C) returns the Rayleigh quotient
%
%                               [-------]         [-------]
%                      --( C )--|   X   |--( C )--|   X   |--( C )--
%                    /          [-------]         [-------]          \
%                    |              |                 |              |
%                  /---\          /---------------------\          /---\
%       theta  =  ( W_L )        (           H           )        ( V_R )
%                  \---/          \---------------------/          \---/
%                    |              |                 |              |
%                    \          [-------]         [-------]          /
%                      --( C )--|   X   |--( C )--[   X   ]--( C )--
%                               [-------]         [-------]
%
%   where WL is the left eigenvector corresponding to the dominant
%   eigenvalue 1 of the left transfer matrix, and VR the right eigenvector
%   corresponding to the dominant eigenvalue 1 of the right transfer
%   matrix.
%
%   iTR1_rq(H,X,C,VR,WL) uses the given matrices VR and WL. Note that VR
%   and WL should be normalized so that trace(WL'*C*VR*C') = 1.
%
%   See also iTR1, iTR1c, iTR1c_rq.

%   Roel Van Beeumen
%   March 18, 2024

if (nargin == 2) || (nargin == 4)
    theta = rq_XX(varargin{:});
elseif (nargin == 3) || (nargin == 5)
    theta = rq_XLCXR(varargin{:});
else
    error('RVB: only 2, 3, 4, or 5 input arguments are allowed!');
end

end


function theta = rq_XX(H,X,V,W)

%% dimensions
[x1,x2,~] = size(X);    assert(x1 == x2);
if nargin == 4
    [v1,v2] = size(V);  assert(v1 == v2); assert(v1 == x1);
    [w1,w2] = size(W);  assert(w1 == w2); assert(w1 == x1);
end
r = x1;

%% default arguments
if nargin < 4
    [eta,V,W] = iTR1_tfmat_domeig(X);
    if abs(eta - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
end

%% convert to matrices
XV = t3_tm(X,V);
XXV = reshape(t4_tt(X,XV),r^2,[]);
WX = t3_mt(W,X);
WXX = reshape(t4_tt(WX,X),r^2,[]);

%% raleigh quotient
theta = sum(sum((WXX'*XXV).*H));

end


function theta = rq_XLCXR(H,X,C,VR,WL)

%% dimensions
[x1,x2,~] = size(X);     assert(x1 == x2);
[c1,c2] = size(C);       assert(c1 == c2); assert(c1 == x1);
if nargin == 5
    [v1,v2] = size(VR);  assert(v1 == v2); assert(v1 == x1);
    [w1,w2] = size(WL);  assert(w1 == w2); assert(w1 == x1);
end
r = x1;

%% left and right cores
XL = t3_mt(C,X);
XR = t3_tm(X,C);

%% default arguments
if nargin < 5
    % left dominant eigenvector
    [etaL,WL] = iTR1_tfmat_domeig(XL,'L');
    if abs(etaL - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
    % right dominant eigenvector
    [etaR,VR] = iTR1_tfmat_domeig(XR,'R');
    if abs(etaR - 1) > 100*eps
        warning('RVB: dominant eigenvalue is not 1!');
    end
    % normalize so that trace(WL'*C*VR*C') = 1
    WL = WL/trace(WL'*C*VR*C');
end

%% convert to matrices
XRVR = t3_tm(XR,VR);
XLCXRVR = reshape(t4_tmt(XL,C,XRVR),r^2,[]);
WLXL = t3_mt(WL,XL);
WLXLCXR = reshape(t4_tmt(WLXL,C,XR),r^2,[]);

%% raleigh quotient
theta = sum(sum((WLXLCXR'*XLCXRVR).*H));

end
