function [theta,X,Y,Sxy,Syx,Theta,Res,Err,idx,wtime] = ...
    iTR2_flexpi(H,tau,A,B,Mab,Mba,maxit,resfreq,maxrank,stagtol,verbose,lamt)
%ITR2_FLEXPI  Flexible power iteration for iTR with 2 cores.
%
%   Inputs:
%     H        Hamiltonian
%     tau      time(s)
%     A        3-tensor
%     B        3-tensor
%     Mab      matrix between A and B
%     Mba      matrix between B and A
%     maxit    maximum number of iterations       [100]
%     resfreq  frequency for computing residual   [1./tau]
%     maxrank  maximum rank of A and B            [max(rank(Mab),rank(Mba))]
%     stagtol  tolerance for residual stagnation  [1e-3]
%     verbose  level of display                   [1]
%     lamt     true eigenvalue (optional)
%
%   Outputs:
%     theta    Rayleigh quotient
%     X        3-tensor
%     Y        3-tensor
%     Sxy      matrix between X and Y
%     Syx      matrix between Y and X
%     Theta    Rayleigh quotients in each iteration
%     Res      residual norm in each iteration
%     Err      truncation error in each iteration
%     idx      indices for each time iteration
%     wtime    wall clock time(s) for each tau
%
%   See also iTR2_pi, iTR2c_flexpi, iTR2c_pi.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 7, maxit = 100; end
if nargin < 8, resfreq = 1./tau; end
if nargin < 9, maxrank = max(rank(Mab),rank(Mba)); end
if nargin < 10, stagtol = 1e-3; end
if nargin < 11, verbose = 1; end
if nargin < 12, lamt = nan; end

%% initialize
X = A;
Y = B;
Sxy = Mab;
Syx = Mba;
k = length(tau);
if isscalar(maxit),   maxit   = maxit*ones(k,1); end
if isscalar(maxrank), maxrank = maxrank*ones(k,1); end
if isscalar(stagtol), stagtol = stagtol*ones(k,1); end
if isscalar(verbose), verbose = verbose*ones(k,1); end
Theta = nan(k+1,1);
Res   = zeros(k+1,1);
Err   = zeros(k+1,1);
idx   = zeros(k+1,1);
wtime = zeros(k,1);

%% flexible power iteration
for i = 1:k
    tstart = tic;
    [theta,X,Y,Sxy,Syx,Theta_tmp,Res_tmp,Err_tmp] = iTR2_pi(H,tau(i),...
        X,Y,Sxy,Syx,maxit(i),resfreq(i),maxrank(i),stagtol(i),verbose(i),lamt);
    wtime(i) = toc(tstart);
    nb = length(Theta_tmp);
    idx(i+1) = idx(i) + nb;
    Theta(idx(i)+1:idx(i+1)) = Theta_tmp;
    Res(idx(i)+1:idx(i+1)) = Res_tmp;
    Err(idx(i)+1:idx(i+1)) = Err_tmp;
end

end
