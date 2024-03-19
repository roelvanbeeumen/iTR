function [X,Y,Sxy,Syx] = iTR2c(A,B,Mab,Mba)
%ITR2C  Canonical form of infinite tensor ring with 2 cores.
%
%                [-------]           [-------]
%     --( Syx )--|   X   |--( Sxy )--|   Y   |--( Syx )--
%   /            [-------]           [-------]            \
%   |                |                   |                |  =  1
%   \            [-------]           [-------]            /
%     --( Syx )--|   X   |--( Sxy )--|   Y   |--( Syx )--
%                [-------]           [-------]
%
%   and
%
%                [-------]           [-------]
%     --( Sxy )--|   Y   |--( Syx )--|   X   |--( Sxy )--
%   /            [-------]           [-------]            \
%   |                |                   |                |  =  1
%   \            [-------]           [-------]            /
%     --( Sxy )--|   Y   |--( Syx )--|   X   |--( Sxy )--
%                [-------]           [-------]
%
%   See also ITR2, ITR1, ITR1C.

%   Reference:
%   R. Orus and G. Vidal. Infinite time-evolving block decimation algorithm
%   beyond unitary evolution. Physical Review B 78, 155117, 2008.
%
%   Roel Van Beeumen
%   March 18, 2024

%% compact
compact = (nargin == 2);

%% dimensions
[a1,a2,n1] = size(A);
[b1,b2,n2] = size(B);
if compact
    assert(a1 == b2);
    assert(a2 == b1);
else
    [r1,r2] = size(Mab);  assert(r1 == a2); assert(r2 == b1);
    [s1,s2] = size(Mba);  assert(s1 == b2); assert(s2 == a1);
end

%% canonical form of big iTR1
if compact
    AB = t3_tt(A,B);
    [Q,Syx] = iTR1c(AB);
else
    AMB = t3_tmt(A,Mab,B);
    [Q,Syx] = iTR1c(AMB,Mba);
end

%% convert back to iTR2
[q1,q2,~] = size(Q);
SQS = t3_mtm(Syx,Q,Syx);
[U,Sxy,V] = t4_svd(reshape(SQS,q1,q2,n1,n2));

%% canonical form
X = t3_mt(inv(Syx),U);
Y = t3_tm(V,inv(Syx));

end
