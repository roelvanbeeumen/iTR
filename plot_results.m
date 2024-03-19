function plot_results(filename,fig)
%PLOT_RESULTS   Plot results
%
%   plot_results(filename) visualizes the results in filename.
%
%   plot_results(filename,fig) visualizes the results in filename in figure fig.
%
%   See also run_ising, run_afm, run_spin.

%   Roel Van Beeumen
%   March 18, 2024

%% default parameter
if nargin < 2, fig = []; end

%% load data
load(filename,'Theta','Res','Err','Idx','tau','lamt','r','docan');

%% plot
if isempty(fig), figure; else, figure(fig); end
idx = 1:Idx(end); idx = idx(isfinite(Theta));
subplot(2,1,1); hold off;
subplot(2,1,2); hold off;
for i = 1:length(Idx)-1
    idx2 = Idx(i)+1:Idx(i+1);
    subplot(2,1,1); semilogx(idx2,Theta(idx2),'*'); hold on;
    subplot(2,1,2); loglog(idx2,abs(Theta(idx2) - lamt),'*'); hold on;
end

subplot(2,1,1);
semilogx(xlim,lamt*[1 1],':k');
lgnd = cell(length(tau)+1,1);
for i = 1:length(tau), lgnd{i} = sprintf('t = %0.e',tau(i)); end
lgnd{end} = '\lambda_t';
legend(lgnd);
% title
if contains(filename,'ising'), prbl = 'Ising Model in a Transverse Field';
elseif contains(filename,'spin'), prbl = 'Heissenberg Spin Model';
elseif contains(filename,'afm'), prbl = 'Heissenberg Isotropic Antiferromagnet';
end
if docan
    title(sprintf('%s   -   iTR2c: r = %i',prbl,r));
else
    title(sprintf('%s   -   iTR2: r = %i',prbl,r));
end
xlabel('iteration'); ylabel('\theta');

subplot(2,1,2);
pres = loglog(idx,Res(idx),'+r');
perr = loglog(Err,':k');
legend([pres,perr],'||residual||','SVD err');
% title
title(filename);
xlabel('iteration'); ylabel('|\theta - \lambda_t|');
tmp = xlim; xlim([0,tmp(2)]);

end
