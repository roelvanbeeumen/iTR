function [theta,X,Y,Sxy,Syx,Theta,Res,Err] = ...
    iTR2c_pi(H,tau,X,Y,Sxy,Syx,maxit,resfreq,maxrank,stagtol,verbose,lamt)
%ITR2C_PI   Power iteration for canonical iTR with 2 cores.
%
%   Inputs:
%     H        Hamiltonian
%     tau      time
%     X        3-tensor
%     Y        3-tensor
%     Sxy      matrix between X and Y
%     Syx      matrix between Y and X
%     maxit    maximum number of iterations       [100]
%     resfreq  frequency for computing residual   [1/tau]
%     maxrank  maximum rank of Sxy and Syx        [max(rank(Sxy),rank(Syx))]
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
%   See also iTR2c_flexpi, iTR2_pi, iTR2_flexpi.

%   Roel Van Beeumen
%   March 18, 2024

%% default arguments
if nargin < 7, maxit = 100; end
if nargin < 8, resfreq = 1/tau; end
if nargin < 9, maxrank = max(rank(Sxy),rank(Syx)); end
if nargin < 10, stagtol = 1e-3; end
if nargin < 11, verbose = 1; end
if nargin < 12, lamt = nan; end

%% dimensions
[x1,x2,~] = size(X);
[y1,y2,~] = size(Y);
[r1,r2] = size(Sxy);  assert(r1 == x2); assert(r2 == y1);
[s1,s2] = size(Syx);  assert(s1 == y2); assert(s2 == x1);

%% initialize
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
        % Rayleigh quotient
        theta = iTR2c_rq(H,X,Y,Sxy,Syx);
        Theta(i) = theta;
        % residual
        res = iTR2c_res(theta,H,X,Y,Sxy,Syx);
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
[U,Sxy,V,e] = t4_svd(HZ,r);
X = t3_mt(inv(Syx),U);
Y = t3_tm(V,inv(Syx));

%% canonical form
[X,Y,Sxy,Syx] = iTR2c(X,Y,Sxy,Syx);

end
