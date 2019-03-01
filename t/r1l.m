%% initialization

rng('default');                               % reset random generator
fun   = @(nodes) prod(bspline_o2(nodes),2);   % example function
d     = 6;                                    % set dimension        
I     = coord_symhc(16,d);                    % choose Index set
z     = [1 33 579 3628 21944 169230 1105193 ...
  77998320 49768670 320144128 2040484044].';  % from Tonis paper
M     = z(d+1);
z     = z(1:d);
x     = zeros(M,d);                           % nodes in space domain
for idx = 1:d
  x(:,idx) = mod(z(idx)*(0:M-1).'/M,1);
end
f     = fun(x);                             % function values
f_hat = bspline_o2_hat(I);                    % original f_hat

f_e   = f+(max(f)-min(f))*0.05*randn(size(f));% noisy function values
W     = 1/M;                                  % weights in space domain
W_hat = max(prod(I,2).^2,1);                  % weights in frequency domain

lambda= 2.^(linspace(-9,0,25));               % possible lambda
err   = 0*lambda;                             % stores l_2-error^2
cv    = 0*lambda;                             % stores LOOCV score


%% main computations

wb = waitbar(0);
for idx = 1:length(lambda) % loop over lambda
  waitbar(idx/length(lambda),wb);
  f_hat_r = compute_f_hat(f_e,I,z,M,lambda(idx),W,W_hat);
  err(idx) = norm(f_hat-f_hat_r);
  f_r = lfft(I,z,M,f_hat_r.');
  h = sum(1./(1+lambda(idx)*W_hat))*W;
  cv(idx) = norm((f_r-f_e)./(1-h))^2;
end
close(wb)

%% computations for plotting

[~,idx_err] = min(err);
[~,idx_cv]  = min(cv);
f_r = F(I,z,M,compute_f_hat(f_e,I,z,M,lambda(idx_cv),W,W_hat));


%% plotting

% plot L_2-error
yyaxis left;
loglog(lambda,err); hold on;
scatter(lambda(idx_err),err(idx_err),'filled'); hold off;
xlabel('\lambda');
ylabel('L_2-error');
% plot cv score
yyaxis right;
loglog(lambda,cv); hold on;
scatter(lambda(idx_cv),cv(idx_cv),'filled'); hold off;
ylabel('cv score');
axis square;


%% helper functions

function y = F(I,z,M,f_hat)
  y = lfft(I,z,M,f_hat.');
end

function f_hat = compute_f_hat(f,I,z,M,lambda,W,W_hat)
  f_hat = alfft(I,z,M,W*f);
  f_hat = f_hat./(1+lambda*W_hat);
end

function f = lfft(I,z,M,f_hat)
  k = mod(I*z,M)+1;
  k = floor(k);
  fhat1 = accumarray(k,f_hat,[M,1],@sum);
  f = M*ifft(fhat1);
end

function f_hat = alfft(I,z,M,f)
  k = mod(I*z,M)+1;
  k = floor(k);
  ghat = fft(f);
  f_hat = ghat(k);
end


function val = bspline_o2( x )
%bspline_o2 quadratic bspline
  x = x - floor(x);

  val = ones(size(x,1),1);

  for t=1:size(x,2)
    ind = find((0<=x(:,t)).*(x(:,t)<1/2));
    if ~isempty(ind)
      val(ind) = val(ind) .* 4.* x(ind,t);
    end

    ind = find((1/2<=x(:,t)).*(x(:,t)<1));
    if ~isempty(ind)
      val(ind) = val(ind) .* 4 .* (1-x(ind,t));
    end

    val = sqrt(3/4) * val;
  end
end

function val = bspline_o2_hat( k )
  %bspline_o2 quadratic bspline
  val = ones(size(k,1),1);

  for t=1:size(k,2)
    ind = find(k(:,t)~=0);
    if ~isempty(ind)
      val(ind) = val(ind) .* bspline_sinc(pi/2*k(ind,t)).^2 .* cos(pi*k(ind,t));
    end

  %   ind = find(k(:,t)==0);
  %   val(ind) = val(ind) .* 1;
    val = sqrt(3/4) * val;
  end
end

function val = bspline_sinc( x )
  if size(x,2) ~= 1
    error('only column vectors are supported');
  end

  val = ones(size(x));
  ind = find(x~=0);
  val(ind) = sin(x(ind))./x(ind);
end
