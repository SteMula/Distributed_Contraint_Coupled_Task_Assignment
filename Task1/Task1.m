%
% DCS MATLAB Project 3, Task 1
% Distributed Task Assignment for Robotic Networks
% Stefano Mulargia, Stefano Bortolotti
% April 2021
%
%%

% init
clear all; 
clc; close all;
addpath .\Functions;

N = 3;              % number of agents
NI = 4;            % size of each agent 
%  !!! Important NI>N  !!!
MAXITERS = 1e4;
S = NI;

%step-size
gamma_const = 3;
gamma_exp = 0.6;


%% Generation matrix + Centralized solution

[c, D, d, H_eq, b_eq, LB, UB] = problem_generator_function(NI,N);
[xopt, fopt] = linprog(c,D,d,H_eq,b_eq,LB,UB);

b_eq = b_eq/N;

%% Network graph
p = 0.2;
[AA_NW, AA] = binomialGraph(p, N, 'doubly');

%% Distributed dual subgradient
options = optimoptions('linprog','Display','none');

% reshape per comodità matrici nei cicli for
c_d = reshape(c,NI,N);
LB_d = reshape(LB,NI,N);
UB_d = reshape(UB,NI,N);
H_eq_d = reshape(H_eq,S,NI,N);
D_d = reshape(nonzeros(D),NI,N)';
D_d1 = reshape(D,N,NI,N);

ZZ = zeros(NI,N,MAXITERS);
ZZ_RA = zeros(NI,N,MAXITERS);

Lambda = zeros(S,N,MAXITERS);
VV = zeros(S,N);

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);
primal_cost_RA = zeros(MAXITERS,1);
dual_cost_RA   = zeros(MAXITERS,1);

consensus_err = zeros(MAXITERS,1);

%%
for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %d\n',tt);
  end
    gamma_t = gamma_const*(1/tt)^gamma_exp; % 1/tt^alpha with alpha in (0.5, 1]

  for ii=1:N
    N_ii = find(AA_NW(:,ii) == 1)';

    VV(:,ii) = AA(ii,ii)*Lambda(:,ii,tt);
    
    for kk = N_ii
      VV(:,ii) = VV(:,ii) + AA(ii,kk)*Lambda(:,kk,tt);
    end
  end

% Primal Update
  for ii=1:N
        [ZZ(:,ii,tt)] = linprog(c_d(:,ii)+(VV(:,ii)'*H_eq_d(:,:,ii))',D_d1(:,:,ii),d,[],[],LB_d(:,ii),UB_d(:,ii),options);
  end

% Running average
  for ii=1:N
      if tt==1
          ZZ_RA(:,ii,tt) = ZZ(:,ii,tt);
      else
          ZZ_RA(:,ii,tt) = (1/tt)*((tt-1)*ZZ_RA(:,ii,tt-1)+ZZ(:,ii,tt));
      end
  end
  
% Dual Update
  for ii=1:N
    grad_ii = H_eq_d(:,:,ii)*ZZ(:,ii,tt)-b_eq;

    Lambda(:,ii,tt+1) = VV(:,ii) + gamma_t*grad_ii;
  end

% Performance check
	Lambda_avg = mean(Lambda(:,:,tt),2);

% Computation of primal cost, dual cost and consensus error
  for ii=1:N
    ff_ii = c_d(:,ii)'*ZZ(:,ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = c_d(:,ii)'*ZZ_RA(:,ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = c_d(:,ii)'*ZZ(:,ii,tt) + Lambda(:,ii,tt)'*(H_eq_d(:,:,ii)*ZZ(:,ii,tt)-b_eq);
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = c_d(:,ii)'*ZZ_RA(:,ii,tt) + Lambda(:,ii,tt)'*(H_eq_d(:,:,ii)*ZZ_RA(:,ii,tt)-b_eq);
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;

    consensus_err(tt) = consensus_err(tt) + norm(Lambda(:,ii,tt) - Lambda_avg);
  end
end

%%
% Last value [for plot]
tt = MAXITERS;
fprintf('Iteration n. %d\n',tt);
Lambda_avg = mean(Lambda(:,tt),2);

for ii=1:N
  N_ii = find(AA_NW(:,ii) == 1)';
  
  VV(:,ii) = AA(ii,ii)*Lambda(:,ii,tt);
  
  for jj = N_ii
    VV(:,ii) = VV(:,ii) + AA(ii,jj)*Lambda(:,jj,tt);
  end
end

for ii=1:N
    [ZZ(:,ii,tt)] = linprog(c_d(:,ii)+(VV(:,ii)'*H_eq_d(:,:,ii))',D_d1(:,:,ii),d,[],[],LB_d(:,ii),UB_d(:,ii),options);

    ZZ_RA(:,ii,tt) = (1/tt)*((tt-1)*ZZ_RA(:,ii,tt-1)+ZZ(:,ii,tt));

    ff_ii = c_d(:,ii)'*ZZ(:,ii,tt);
    primal_cost(tt) = primal_cost(tt) + ff_ii;
    
    ff_ii = c_d(:,ii)'*ZZ_RA(:,ii,tt);
    primal_cost_RA(tt) = primal_cost_RA(tt) + ff_ii;
    
    qq_ii = c_d(:,ii)'*ZZ(:,ii,tt) + Lambda(:,ii,tt)'*(H_eq_d(:,:,ii)*ZZ(:,ii,tt)-b_eq);
    dual_cost(tt) = dual_cost(tt) + qq_ii;
    
    qq_ii = c_d(:,ii)'*ZZ_RA(:,ii,tt) + Lambda(:,ii,tt)'*(H_eq_d(:,:,ii)*ZZ_RA(:,ii,tt)-b_eq);
    dual_cost_RA(tt) = dual_cost_RA(tt) + qq_ii;
end

%% dual cost
figure
  semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(dual_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('dual cost error')
  legend('dual cost','dual cost with running avg')

%% primal cost
figure
  semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(primal_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('primal cost error')
  legend('primal cost','primal cost with running avg')

%% primal vs dual
% figure
%   semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
%   hold on, grid on, zoom on
%   semilogy(1:MAXITERS,abs(primal_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
%   semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
%   semilogy(1:MAXITERS,abs(dual_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
%   xlabel('t')
%   ylabel('cost error')
%   legend('primal cost','primal cost with running avg','dual cost','dual cost with running avg')

%% lambda
temp=reshape(Lambda,NI*N,tt);
figure
  plot(1:MAXITERS,temp, 'LineWidth',1.5);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('\lambda_i^t')

%% plot consensus error
figure
  semilogy(1:MAXITERS,consensus_err(1:MAXITERS), 'LineWidth',2);
  hold on, grid on, zoom on
  refresh
  xlabel('t')
  ylabel('consensus error')
