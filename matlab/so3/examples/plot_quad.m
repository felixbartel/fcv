%
% this script calculates and plots the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

rng('default');                           % reset random generatormber of nodes
N         = 23;                           % bandwidth
N = floor(N/2);
odf = SantaFe;
fhat = odf.calcFourier;
fhat = fhat(1:(N+1)*(4*(N+1)*(N+1)-1)/3);
odf = FourierODF(fhat,crystalSymmetry('m-3m'),specimenSymmetry('222'));

fun       = @(alpha,beta,gamma) odf.eval(...
  orientation.byEuler(alpha,beta,gamma));    % example function

A = load('examples/N23_M5880_IcoC7.dat');
nodes = orientation.byMatrix(reshape(A',3,3,[]));
nodes = nodes(:);
[alpha,beta,gamma] = nodes.Euler;

M         = length(nodes);                % number of nodes
W = 8*pi^2/M;
f         = fun(alpha,beta,gamma);        % function values

f_e       = f+0.1*randn(size(f))*(max(f)-min(f)); % noisy function values
s         = 3;                            % weights in requency domain

lambda    = 2.^(linspace(-25,-18,25));    % possible lambda
err       = 0*lambda;                     % stores L_2-error
cv        = 0*lambda;                     % stores cv score


%% main computations

fcv = FCV_quad([alpha,beta,gamma],f_e,W,N,s);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  [cv(idx),~,fhat_r] = fcv.compute(lambda(idx));
  err(idx) = norm(fhat-fhat_r);
end
close(wb);


%% computations for plotting

% get optimal lambda
[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

[~,~,fhat_r] = fcv.compute(lambda(idx_cv));
odf_r = FourierODF(fhat_r,crystalSymmetry('m-3m'),specimenSymmetry('222'));

[~,~,fhat_r] = fcv.compute(0);
odf_noisy = FourierODF(fhat_r,crystalSymmetry('m-3m'),specimenSymmetry('222'));

%% plotting

figure(1);
% plot L_2-error
yyaxis left;
loglog(lambda,err); hold on;
scatter(lambda(idx_err),err(idx_err),40,'filled'); hold off;
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
p = loglog(lambda,cv); hold on;
scatter(lambda(idx_cv),cv(idx_cv),40,'filled'); hold off;
ylabel('cv score');
axis square;

figure(2);
plot(odf,'sigma');
mtexColorbar;
mtexColorMap(jet);

figure(3);
plot(odf_r,'sigma');
mtexColorbar;
mtexColorMap(jet);

figure(4);
plot(odf_noisy,'sigma');
mtexColorbar;
mtexColorMap(jet);
