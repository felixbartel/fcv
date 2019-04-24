%
% this script calculates and plots the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

rng('default');                           % reset random generatormber of nodes
N         = 23;                           % bandwidth
N = floor(N/2);
mtexdata 'ol';
odf = calcODF(ebsd('ol').orientations);
fhat = odf.calcFourier;                   % oroginal fhat
fhat = [fhat;zeros((N+1)*(4*(N+1)*(N+1)-1)/3-length(fhat),1)];
fhat = fhat(1:(N+1)*(4*(N+1)*(N+1)-1)/3);
odf = FourierODF(fhat,odf.CS,odf.SS);

A = load('examples/N23_M5880_IcoC7.dat');
nodes = orientation.byMatrix(reshape(A',3,3,[]));
nodes = nodes(:);
angles = Euler(nodes,'nfft');

M         = length(nodes);                % number of nodes
W = 8*pi^2/M;
f         = odf.eval(nodes);        % function values

f_e       = f+0.1*randn(size(f))*(max(f)-min(f)); % noisy function values
s         = 3;                            % weights in requency domain

lambda    = 2.^(linspace(-18,-15,25));    % possible lambda
err       = 0*lambda;                     % stores L_2-error
cv        = 0*lambda;                     % stores cv score


%% main computations

fcv = FCV_quad(angles,f_e,W,N,s);

n = (0:N)';
n = repelem(n,(1:2:(2*N+1)).^2);

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  [cv(idx),~,fhat_r] = fcv.compute(lambda(idx));
  err(idx) = norm((fhat-fhat_r).*(2*n+1)/(8*pi^2));
end
close(wb);


%% computations for plotting

% get optimal lambda
[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);

[~,~,fhat_r] = fcv.compute(lambda(idx_cv));
odf_r = FourierODF(fhat_r,odf.CS,odf.SS);

nodes.CS = odf.CS;
nodes.SS = odf.SS;


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
plotSection(nodes,f,'sigma','all','MarkerSize',1);
mtexColorbar;
mtexColorMap('WhiteJet');
c = caxis();

figure(3);
plotSection(nodes,f_e,'sigma','all','MarkerSize',1);
mtexColorbar;
mtexColorMap('WhiteJet');
setColorRange(c);

figure(4);
plot(odf_r,'sigma');
mtexColorbar;
mtexColorMap('WhiteJet');
setColorRange(c);

figure(5);
plot(odf,'sigma');
mtexColorbar;
mtexColorMap('WhiteJet');
setColorRange(c);

