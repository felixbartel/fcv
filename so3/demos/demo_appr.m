%
% this script calculates and plots the cv score for an quadrature
% example on the two-dimensional sphere
%
%% initialization

rng('default');                           % reset random generatormber of nodes

mtexdata 'ol';
odf = calcODF(ebsd('ol').orientations);
N         = odf.bandwidth-1;              % bandwidth
fhat = odf.calcFourier;                   % oroginal fhat
fhat = [fhat;zeros((N+1)*(4*(N+1)*(N+1)-1)/3-length(fhat),1)];
fhat = fhat(1:(N+1)*(4*(N+1)*(N+1)-1)/3);
odf = FourierODF(fhat,odf.CS,odf.SS);

M         = 2*(N+1)*(4*(N+1)^2-1)/3;      % number of nodes
nodes = orientation.rand(M);
angles = Euler(nodes,'nfft');
W = VoronoiArea(nodes);

f         = odf.eval(nodes);              % function values

f_e       = f+0.05*randn(size(f))*(max(f)-min(f)); % noisy function values
s         = 3;                            % weights in requency domain

lambda    = 2.^(linspace(-22,-20,25));     % possible lambda
err       = 0*lambda;                     % stores L_2-error
ocv       = 0*lambda;                     % stores ocv score
gcv       = 0*lambda;                     % stores gcv score
ocv_exact = 0*lambda;                     % stores exact ocv score
gcv_exact = 0*lambda;                     % stores exact gcv score


%% main computations

fcv = FCV_appr(angles,f_e,W,N,s);

n = (0:N)';
n = repelem(n,(1:2:(2*N+1)).^2);
tic
wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  res = fcv.compute(lambda(idx));
  ocv(idx) = res.ocv;
  gcv(idx) = res.gcv;
%  res = fcv.compute(lambda(idx),"exact");
%  ocv_exact(idx) = res.ocv;
%  gcv_exact(idx) = res.gcv;
  err(idx) = norm((fhat-res.fhat_r).*(2*n+1)/(8*pi^2));
end
close(wb);
toc

%% computations for plotting

% get optimal lambda
[~,idx_err] = min(err);
[~,idx_gcv_exact] = min(gcv_exact);
[~,idx_ocv_exact] = min(ocv_exact);
[~,idx_gcv] = min(gcv);
[~,idx_ocv] = min(ocv);

res = fcv.compute(lambda(idx_gcv));
odf_r = FourierODF(res.fhat_r,odf.CS,odf.SS);

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
p = loglog(lambda,[ocv_exact; gcv_exact; ocv; gcv]); hold on;
scatter(lambda([idx_ocv_exact idx_gcv_exact idx_ocv idx_gcv]),...
  [ocv_exact(idx_ocv_exact) gcv_exact(idx_gcv_exact) ocv(idx_ocv) gcv(idx_gcv)],40,'filled');
ylim([min(real([ocv_exact gcv_exact gcv])) max(real([ocv_exact gcv_exact gcv]))]); hold off;
legend(p,'ocv','gcv','appr ocv','appr gcv');
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

