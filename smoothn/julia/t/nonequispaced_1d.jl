#
# this script calculates and plots the cv score for an nonequispaced
# example on the one-dimensional torus
#
push!(LOAD_PATH, string(@__DIR__, "/../../../../nfft/julia/nfft"))
using NFFT
using IterativeSolvers
using LinearMaps
using LinearAlgebra

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

function gcv(plan, lambda, W, W_hat, f)
  f_r = F(plan, compute_f_hat(plan, lambda, W, W_hat, f));
  h = sum(1 ./(1 .+lambda*W_hat))*W;
  h = sum(h)/length(h);
  return norm((f_r-f)./(1 .-h))^2;
end

function smoothn(nodes, f_e, N_iter = 20)
  nodes     = sort(nodes);
  M         = length(nodes);            # number of nodes
  N         = floor(Int, M/2);          # bandwidth
  s         = 3;                        # weights in frequency domain
  W_hat     = [ 1+abs(n).^s for n in -N/2:N/2-1 ];

## calculate Voronoi area (weights in space domain)
  W = zeros(size(nodes));
  W[1] = nodes[2]-nodes[end]+1;
  W[2:end-1] = nodes[3:end]-nodes[1:end-2];
  W[end] = nodes[1]+1-nodes[end-1];
  W = W/2;

## main computations
  plan = Plan((N, ), M);
  plan.x = nodes;

  fun(x) = log(gcv(plan, exp(x), W, W_hat, f_e));

  deltax = 0.01;

  x_l = -25;
  x_r = 0;
  x_m = (x_l+x_r)/2;
  D_l = (fun(x_l+deltax)-fun(x_l-deltax))/2deltax;
  D_r = (fun(x_r+deltax)-fun(x_r-deltax))/2deltax;
  D_m = (fun(x_m+deltax)-fun(x_m-deltax))/2deltax;

  for i = 1:N_iter
    if D_l*D_m > 0
      x_l = x_m;
      D_l = D_m;
      x_m = (x_l+x_r)/2;
      D_m = (fun(x_m+deltax)-fun(x_m-deltax))/2deltax;
    elseif D_r*D_m > 0
      x_r = x_m;
      D_r = D_m;
      x_m = (x_l+x_r)/2;
      D_m = (fun(x_m+deltax)-fun(x_m-deltax))/2deltax;
    end
    println([ x_m D_m ]);
    if abs(D_m) < 1e-10
      break;
    end
  end
  return compute_f_hat(plan, exp(x_m), W, W_hat, f);
end


using Plots

nodes = rand(Float64, 1000);
nodes = sort(nodes);
f = [ sin(pi*4*x)/2+0.2*sin(pi*32*x)+0.5+0.05*randn() for x in nodes ];
f_hat_r = smoothn(nodes, f);

M_plot = floor(Int, 1e4);
nodes_plot = collect((0:M_plot-1)/M_plot);
plan = Plan((length(f_hat_r), ), M_plot);
plan.x = nodes_plot;
f_r = F(plan, f_hat_r);


scatter(nodes, real(f), markersize = 1, markercolor = :black)
plot!(nodes_plot, real(f_r), linewidth = 2, color = :blue)
