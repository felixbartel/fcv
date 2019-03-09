#
# this script calculates and plots the cv score for an nonequispaced
# example on the one-dimensional torus
#
## helper functions

function F(plan, f_hat)
  plan.fhat = f_hat;
  NFFT.trafo(plan)
  return copy(plan.f);
end

function A_noconj(plan, lambda, W, W_hat, f_hat)
  plan.fhat = f_hat;
  NFFT.trafo(plan);
  y = (W.^0.5).*plan.f;
  return append!(y, ((lambda*W_hat).^0.5).*f_hat);
end

function A_conj(plan, lambda, W, W_hat, f)
  y = (W.^0.5).*f[1:length(W)];
  plan.f = y;
  NFFT.adjoint(plan);
  return plan.fhat+((lambda*W_hat).^0.5).*f[length(W)+1:end];
end

function compute_f_hat(plan, lambda, W, W_hat, f)
  M = length(W);
  N = length(W_hat);
  f = (W.^0.5).*f;
  append!(f, zeros(length(W_hat), 1));
  A = LinearMap(x->A_noconj(plan, lambda, W, W_hat, x), x->A_conj(plan, lambda, W, W_hat, x), M+N, N);
  return lsqr(A, complex(f), maxiter = 20);
end


## initialization

push!(LOAD_PATH, string(@__DIR__, "/../../../../nfft/julia/nfft"))
using FFTW
using NFFT
using IterativeSolvers
using LinearMaps
using LinearAlgebra
using PyPlot

fun(x)    = (3*(1-x).^2 .*exp(-((6x-3).^2)-1)-10*((6x-3)/5-(6x-3).^3).*exp(-(6x-3).^2)-1/3*exp(-((6x-3)+1).^2)+4)/7; # function
N         = 64;                       # bandwidth
nodes     = rand(Float64, 2N).^2;      # nodes
nodes     = sort(nodes);
M         = length(nodes);            # number of nodes
f         = [fun(n) for n in nodes];  # function values

f_e       = f+0.05*randn(size(f));    # noisy function values
s         = 3;                        # weights in frequency domain
W_hat     = [ abs(n).^s+1 for n in -N/2:N/2-1 ];

lambda    = 2.0 .^(range(-30, length=50, stop=0));
err       = zeros(length(lambda));    # stores L_2-error
ocv       = zeros(length(lambda));    # stores ocv score
gcv       = zeros(length(lambda));    # stores gcv score
ocv_appr  = zeros(length(lambda));    # stores approximated ocv score
gcv_appr  = zeros(length(lambda));    # stores approximated gcv score


## original f_hat

M2 = 2^11;                            # sufficiently large bandwidth
nodes2 = (0:M2-1)/M2;                 # nodes in space domain
fe = [fun(n) for n in nodes2];        # function values
f_hat = ifft(fe);
f_hat = [ f_hat[floor(Int, idx)] for idx = [(length(f_hat)-N/2+1):length(f_hat); 1:N/2] ];


## calculate Voronoi area (weights in space domain)

W = zeros(size(nodes));
W[1] = nodes[2]-nodes[end]+1;
W[2:end-1] = nodes[3:end]-nodes[1:end-2];
W[end] = nodes[1]+1-nodes[end-1];
W = W/2;


## main computations

plan = Plan((N, ), M);
plan.x = nodes;
A(lambda) = LinearMap(x->A_noconj(plan, x, lambda, W, W_hat), x->A_conj(plan, x, lambda, W, W_hat), M+N, N);

for idx = 1:length(lambda) # loop over lambda
  print(string("\r", floor(Int, 100idx/length(lambda)), " %"));
  
  f_hat_r = compute_f_hat(plan, lambda[idx], W, W_hat, f_e);
  f_r = F(plan, f_hat_r);
  err[idx] = norm(f_hat-f_hat_r);
  
# approximated cv score
  h = sum(1 ./(1 .+lambda[idx]*W_hat))*W;
  ocv_appr[idx] = norm((f_r-f_e)./(1 .-h))^2;
  h_mean = sum(h)/length(h);
  gcv_appr[idx] = norm((f_r-f_e)/(1-h_mean))^2;
  
# # exact cv score
#   h = zeros(Complex{Float64}, M);
#   for l = 1:M
#     tmp = zeros(M); tmp[l] = 1;
#     tmp = F(plan, compute_f_hat(plan, lambda[idx], W, W_hat, tmp));
#     h[l] = real(tmp[l]);
#   end
#   ocv[idx] = norm((f_r-f_e)./(1 .-h))^2;  
#   h_mean = sum(h)/length(h);
#   gcv[idx] = norm((f_r-f_e)/(1-h_mean))^2;
end


## computations for plotting

idx_err       = argmin(err);
idx_gcv       = argmin(gcv);
idx_ocv       = argmin(ocv);
idx_gcv_appr  = argmin(gcv_appr);
idx_ocv_appr  = argmin(ocv_appr);

f_hat_r = compute_f_hat(plan, lambda[idx_ocv_appr], W, W_hat, f_e);

res = 480;
plotnodes = collect(range(0, length=res, stop=1));

plan = Plan((N, ), res);
plan.x = plotnodes;
plotf_r = F(plan, f_hat_r);


## plotting

clf();
# plot reconstruction
subplot(121);
scatter(nodes, real(f_e));
plot(plotnodes, real(plotf_r), color = "orange");


# plot L_2-error
ax = subplot(122);
ax2 = ax.twinx()
ax.plot(lambda, err);
ax.scatter([lambda[idx_err]], [err[idx_err]]);

# plot cv score
ax2.plot(lambda, ocv, color = "orange", linestyle = "-");
ax2.plot(lambda, gcv, color = "orange", linestyle = "-.");
ax2.plot(lambda, ocv_appr, color = "orange", linestyle = "--");
ax2.plot(lambda, gcv_appr, color = "orange", linestyle = ":");
ax2.scatter(lambda[[idx_ocv idx_gcv idx_ocv_appr idx_gcv_appr]],
  [ocv[idx_ocv] gcv[idx_gcv] ocv_appr[idx_ocv_appr] gcv_appr[idx_gcv_appr]],
  color="orange");

ax.set_yscale("log");
ax.set_xscale("log");
ax.set_xlabel("lambda");
ax.set_ylabel("L2-error");

ax2.set_yscale("log");
ax2.set_xscale("log");
ax2.set_ylabel("cv-score");
ax2.set_ylim(min([ocv; gcv; gcv_appr]...), max([ocv; gcv; gcv_appr]...));

legend(["ocv", "gcv", "appr ocv", "appr gcv"]);

show()
