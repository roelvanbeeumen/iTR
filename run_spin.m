function filename = run_spin(varargin)
%RUN_SPIN   Heissenberg Spin Model
%
%   run_spin(r,tau,maxit,resfreq,stagtol,verbose,docan)
%
%       r         iTR rank                           [30]
%       tau       time(s)                            [1e-1,1e-2,1e-3]
%       maxit     maximum number of iterations       [1e6]
%       resfreq   frequency for computing residual   [1./tau]
%       stagtol   tolerance for residual stagnation  [1e-3]
%       verbose   level of display                   [1]
%       docan     iTR or iTRc                        [false]
%
%   run_spin(filename,tau,maxit,resfreq,stagtol,verbose)
%
%       filename  checkpoint
%       tau       extra time(s)
%       maxit     maximum number of iterations       [1e6]
%       resfreq   frequency for computing residual   [1./tau]
%       stagtol   tolerance for residual stagnation  [1e-3]
%       verbose   level of display                   [1]

%   Roel Van Beeumen
%   March 18, 2024

%% parameters
newrun = nargin < 1 || isnumeric(varargin{1});
if newrun
    % new run
    if nargin < 1,       r = 30;               else,       r = varargin{1}; end
    if nargin < 2,     tau = [1e-1,1e-2,1e-3]; else,     tau = varargin{2}; end
    if nargin < 3,   maxit = 1e6;              else,   maxit = varargin{3}; end
    if nargin < 4, resfreq = 1./tau;           else, resfreq = varargin{4}; end
    if nargin < 5, stagtol = 1e-3;             else, stagtol = varargin{5}; end
    if nargin < 6, verbose = 1;                else, verbose = varargin{6}; end
    if nargin < 7,   docan = false;            else,   docan = varargin{7}; end
    % initialize
    rng('default');
    A = iTR1(r);
    B = A;
    [A,B,Mab,Mba] = iTR2_normalize(A,B,eye(r),eye(r));
    [X,Y,Sxy,Syx] = iTR2c(A,B,Mab,Mba);
else
    % checkpoint restart
    checkpoint = load(varargin{1},'r','tau','stagtol','docan',...
        'X','Y','Sxy','Syx','Theta','Res','Err','Idx','wtime');
                         r = checkpoint.r;
                       tau = varargin{2};
    if nargin < 3,   maxit = 1e6;              else,   maxit = varargin{3}; end
    if nargin < 4, resfreq = 1./tau;           else, resfreq = varargin{4}; end
    if nargin < 5, stagtol = 1e-3;             else, stagtol = varargin{5}; end
    if nargin < 6, verbose = 1;                else, verbose = varargin{6}; end
                     docan = checkpoint.docan;
    % initialize
    X = checkpoint.X;
    Y = checkpoint.Y;
    Sxy = checkpoint.Sxy;
    Syx = checkpoint.Syx;
end

%% Pauli matrices
sigmax = [0 1; 1 0];
sigmay = [0 -1i; 1i 0];
sigmaz = [1 0; 0 -1];

%% Heissenberg spin model
H = kron(sigmax/2,sigmax/2) + kron(sigmay/2,sigmay/2) + kron(sigmaz/2,sigmaz/2);
lamt = -log(2) + 0.25;

%% flexible power iteration
if docan
    [theta,X,Y,Sxy,Syx,Theta,Res,Err,Idx,wtime] = ...
        iTR2c_flexpi(H,tau,X,Y,Sxy,Syx,maxit,resfreq,r,stagtol,verbose,lamt);
else
    [theta,X,Y,Sxy,Syx,Theta,Res,Err,Idx,wtime] = ...
        iTR2_flexpi(H,tau,X,Y,Sxy,Syx,maxit,resfreq,r,stagtol,verbose,lamt);
end

%% update outputs
if ~newrun
    tau     = [checkpoint.tau    ,tau];
    stagtol = [checkpoint.stagtol,stagtol];
    Theta   = [checkpoint.Theta  ;Theta];
    Res     = [checkpoint.Res    ;Res];
    Err     = [checkpoint.Err    ;Err];
    Idx     = [checkpoint.Idx    ;Idx(2:end) + checkpoint.Idx(end)];
    wtime   = [checkpoint.wtime  ;wtime];
end

%% save
filename = sprintf('spin-r%i-%6s.mat',r,datestr(now,'yymmdd-HHMMSS'));
save(filename,'varargin','r','tau','stagtol','docan',...
    'H','lamt','theta','X','Y','Sxy','Syx','Theta','Res','Err','Idx','wtime');

%% plot
plot_results(filename);

end
