function theta = iTR1c_rq(varargin)
%ITR1C_RQ   Rayleigh quotient of canonical infinite tensor ring with 1 core.
%
%   iTR1c_rq(H,Q,S) uses the canonical form:
%
%                            [-------]         [-------]
%                   --( S )--|   Q   |--( S )--|   Q   |--( S )--
%                 /          [-------]         [-------]          \
%                 |              |                 |              |
%                 |            /---------------------\            |
%       theta  =  |           (           H           )           |
%                 |            \---------------------/            |
%                 |              |                 |              |
%                 \          [-------]         [-------]          /
%                   --( S )--|   Q   |--( S )--|   Q   |--( S )--
%                            [-------]         [-------]
%
%   iTR1c_rq(H,QL,C,QR) uses the mixed canonical form:
%
%                       [-------]         [-------]
%                   ----|  Q_L  |--( C )--|  Q_R  |----
%                 /     [-------]         [-------]     \
%                 |         |                 |         |
%                 |       /---------------------\       |
%       theta  =  |      (           H           )      |
%                 |       \---------------------/       |
%                 |         |                 |         |
%                 \     [-------]         [-------]     /
%                   ----|  Q_L  |--( C )--|  Q_R  |----
%                       [-------]         [-------]
%
%   See also iTR1, iTR1c, iTR1_rq.

%   Roel Van Beeumen
%   March 18, 2024

if nargin == 3
    theta = rq_SQSQS(varargin{:});
elseif nargin == 4
    theta = rq_QLCQR(varargin{:});
else
    error('RVB: only 3 or 4 input arguments are allowed!');
end

end


function theta = rq_SQSQS(H,Q,S)

%% dimensions
[q1,q2,n] = size(Q);
[s1,s2] = size(S);
assert(q1 == s2); assert(q2 == s1);

%% convert to matrix
Z = t4_mtmtm(S,Q,S,Q,S);
Z = reshape(Z,s1*s2,n*n);

%% rayleigh quotient
ZtZ = Z'*Z;
theta = sum(sum(ZtZ.*H));

end


function theta = rq_QLCQR(H,QL,C,QR)

%% dimensions
[c1,c2] = size(C);
[~,l2,n1] = size(QL);  assert(l2 == c1);
[r1,~,n2] = size(QR);  assert(r1 == c2); assert(n1 == n2);

%% convert to matrix
Z = t4_tmt(QL,C,QR);
Z = reshape(Z,[],n1*n2);

%% rayleigh quotient
ZtZ = Z'*Z;
theta = sum(sum(ZtZ.*H));

end
