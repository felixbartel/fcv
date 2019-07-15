M = 8;
M_plot = 2^9;
N = M;

nodes = (0:M-1)'/M;
%nodes = mod(nodes+0.02*randn(size(nodes)),1);

f = double((0:M-1) == M/2)';

What = 1+abs(-N/2:N/2-1).^20;

What = What(:);
fcv = FCV_appr(nodes,f,[],What);

clf; hold on;
scatter(nodes,f,100,'k');
scatter(nodes(M/2+1),f(M/2+1),100,'filled','k');
colormap autumn;

for lambda = 2.^(linspace(-70,10,100))
  
  s = fcv.compute(lambda);
  fhat_r = s.fhat_r;
  hexact = s.f_r;
  hexact = hexact(M/2+1);
  
  f_r = fftd(...
    [fhat_r(end/2+1:end); ...
    zeros((M_plot-length(fhat_r)),1);...
    fhat_r(1:end/2)],1);
  
  h = fcv.W(M/2+1)*sum(1./(1+lambda*What));
  
  plot((0:M_plot-1)/M_plot,real(f_r),'color',0.5*ones(1,3));
%  scatter(nodes(M/2+1),h,50,norm(h-hexact));
  
  drawnow();
end

scatter(nodes,f,100,'k');
scatter(nodes(M/2+1),f(M/2+1),100,'filled','k')

hold off;
