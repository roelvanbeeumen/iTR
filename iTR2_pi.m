function [theta,X,Y,Sxy,Syx,Theta,Res,Err] = ...
    iTR2_pi(H,tau,A,B,Mab,Mba,maxit,resfreq,maxrank,stagtol,verbose,lamt)
%ITR2_PI  Power iteration for iTR with 2 cores.
%
%   Inputs:
%     H        Hamiltonian
%     tau      time
%     A        3-tensor
%     B        3-tensor
%     Mab      matrix between A and B
%     Mba      matrix between B and A
%     maxit    maximum number of iterations       [100]
%     resfreq  frequency for computing residual   [1/tau]
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
%     Theta    Rayleigh quotient in each iteration
%     Res      residual norm in each iteration
%     Err      truncation error in each iteration
%
%   See also iTR2_flexpi, iTR2c_pi, iTR2c_flexpi.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 7, maxit = 100; end
if nargin < 8, resfreq = 1/tau; end
if nargin < 9, maxrank = max(rank(Mab),rank(Mba)); end
if nargin < 10, stagtol = 1e-3; end
if nargin < 11, verbose = 1; end
if nargin < 12, lamt = nan; end

%% dimensions
[a1,a2,~] = size(A);
[b1,b2,~] = size(B);
[r1,r2] = size(Mab);  assert(r1 == a2); assert(r2 == b1);
[s1,s2] = size(Mba);  assert(s1 == b2); assert(s2 == a1);

%% initialize
X = A;
Y = B;
Sxy = Mab;
Syx = Mba;
Theta = nan(maxit,1);
Res   = nan(maxit,1);
Err   = nan(maxit,1);
if verbose
    fprintf('\ntau = %0.e\n\n',tau);
    fprintf('       it         theta            res        SVD err');
    frmt = '%9i  %18.10e  %11.4e  %11.4e';
    if isfinite(lamt)
        fprintf('        diff');
        frmt = [frmt,'  %11.4e'];
    end
    fprintf('\n\n');
end

%% expm(-tau*H)
expH = expm(-tau*H);

%% iteration
thetaopt = nan; resopt = inf; stag = 0;
Xopt = []; Yopt = []; Sxyopt = []; Syxopt = [];
for i = 1:maxit
    
    % X -- Sxy -- Y
    [X,Y,Sxy,Syx,e1] = oper(expH,X,Y,Sxy,Syx,maxrank);
    
    % Y -- Syx -- X
    [Y,X,Syx,Sxy,e2] = oper(expH,Y,X,Syx,Sxy,maxrank);
    
    % SVD truncation error
    Err(i) = max(e1,e2);
    
    if (mod(i,resfreq) == 0) || (i == maxit)
        % canonical decomposition
        [Q,U,S,T] = iTR2c(X,Y,Sxy,Syx);
        % Rayleigh quotient
        theta = iTR2c_rq(H,Q,U,S,T);
        Theta(i) = theta;
        % residual
        res = iTR2c_res(theta,H,Q,U,S,T);
        Res(i) = norm(res(:));
        % output
        if verbose && isfinite(lamt)
            fprintf(frmt,i,theta,Res(i),Err(i),abs(theta - lamt));
        elseif verbose
            fprintf(frmt,i,theta,Res(i),Err(i));
        end
        % check stagnation of residual
        if stagtol ~= 0
            isstag = (resopt - Res(i))/10^(floor(log10(Res(i)))) < stagtol;
            if (Res(i) > resopt) || isstag
                stag = stag + 1;
                if verbose && (Res(i) > resopt), fprintf('   --> no decrease');
                elseif verbose,               fprintf('   --> stagnation'); end
                if stag > 2
                    Theta = Theta(1:i);
                    Res   = Res(1:i);
                    Err   = Err(1:i);
                    if verbose, fprintf('\n'); end
                    break
                end
            else
                stag = 0;
            end
        end
        % update
        if Res(i) < resopt
            thetaopt = theta; resopt = Res(i);
            Xopt = X; Yopt = Y; Sxyopt = Sxy; Syxopt = Syx;
        end
        if verbose, fprintf('\n'); end
    end
    
end

%% return optimal
theta = thetaopt; X = Xopt; Y = Yopt; Sxy = Sxyopt; Syx = Syxopt;

%% output
if verbose && isfinite(lamt)
    fprintf('\n  ==>  |theta - lamt| = %d\n',abs(theta - lamt));
end

end


function [X,Y,Sxy,Syx,e] = oper(expH,X,Y,Sxy,Syx,r)

%% apply exp(-tau*H)
Z = t4_mtmtm(Syx,X,Sxy,Y,Syx);
HZ = t4_oper(expH,Z);

%% split
[X,Sxy,Y,e] = t4_svd(HZ,r);

%% scale
Sxy = Sxy/norm(Sxy,'fro');

%% update
X = t3_mt(inv(Syx),X);
Y = t3_tm(Y,inv(Syx));

end
