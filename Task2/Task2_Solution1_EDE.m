%
% DCS MATLAB Project 3, Task 2, Solution 1 (Everyone can Do Everything)
% Distributed Task Assignment for Robotic Networks
% Stefano Mulargia, Stefano Bortolotti
% April 2021
%
%%
% In this code we developed a task assignment in which every agent can
% reach every task.
%
%%
clear all; 
clc; close all;
addpath .\Functions;
rng(1)
N = 3;          % is the number of Agents
NI = 3;         % size of each agent 
%  !!! Important NI>N  !!!
MAXITERS = 1e4;

Threshold = 0.9;

% Note that if you decreasse the number of iteration you have to decrease
% also the threshold. Unfortunately there is not a formula to compute it
% but you have to make a trade off between number of iteration and the
% threshold.

options = optimoptions('linprog','Display','none');

%% Generation matrix + Centralized solution
b_eq = ones(N,1);
b_eq = b_eq/N;
LB = zeros(N*NI,1);
UB = ones(N*NI,1);
H_eq = zeros(NI,NI*N);

index = 0;
for i=1:N
    H_eq(:,index+1 :index+N) = eye(N);
    index=index +N;
end

G = ones(1,N);
g = 1;

c = zeros(NI*N,1);

%% Network graph
p = 0.5;
[AA_NW, AA] = binomialGraph(p, N, 'doubly');

%% Initialization
dd=3; % dimension (2 if you want 2D agents, 3 if you want 3D agents)

distance = zeros(N,N);

task = rand(dd,N);

xx_init = rand(dd,N);
for i=1:N
    for j=1:N
    distance(i,j) = norm(xx_init(:,i)-task(:,j));
    end
end

% stack equality constraints together
Equal = vertcat(H_eq,blkdiag(G,G,G));
b_equal = vertcat(b_eq*N,ones(N,1));

% copy NxN distance matrix into 3Nx1 cost column vector
index = 0;
for i=1:N
    c(index+1:index+N,1) = distance(:,i);
    index = index + N;
end

%% Centralized Solution
clc
[xopt, fopt] = linprog(c',[],[],Equal,b_equal,LB,UB,options)

%% Distributed Solution

c_d = distance;

LB_d = reshape(LB,NI,N);
UB_d = reshape(UB,NI,N);

H_eq_d = reshape(H_eq,NI,NI,N);
b_eq_d = reshape(b_eq,NI,1);
G_eq_d = blkdiag(G,G,G);
g_d = vertcat(eye(N),eye(N));

Task_assigned = zeros(N,1);
ZZ = zeros(NI,N,MAXITERS);
ZZ_RA = zeros(NI,N,MAXITERS);

Lambda = zeros(NI,N,MAXITERS);
VV = zeros(NI,N);

primal_cost = zeros(MAXITERS,1);
dual_cost   = zeros(MAXITERS,1);
primal_cost_RA = zeros(MAXITERS,1);
dual_cost_RA   = zeros(MAXITERS,1);

gamma_const = 3;
gamma_exp = 0.6;

consensus_err = zeros(MAXITERS,1);
xx_actual_position = xx_init;

%%
for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %g\n',tt);
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
        [ZZ(:,ii,tt)]  = linprog(c_d(:,ii)+H_eq_d(:,:,ii)'*VV(:,ii),[],[],...
            G,g,LB_d(:,ii),UB_d(:,ii),options);
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
	MU_avg = mean(Lambda(:,:,tt),2);

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

    consensus_err(tt) = consensus_err(tt) + norm(Lambda(:,ii,tt) - MU_avg);
  end
end


%%
% Last value [for plot]
tt = MAXITERS;
fprintf('Iteration n. %g\n',tt);
MU_avg = mean(Lambda(:,tt),2);

for ii=1:N
  N_ii = find(AA_NW(:,ii) == 1)';
  
  VV(:,ii) = AA(ii,ii)*Lambda(:,ii,tt);
  
  for jj = N_ii
    VV(:,ii) = VV(:,ii) + AA(ii,jj)*Lambda(:,jj,tt);
  end
end

for ii=1:N

    [ZZ(:,ii,tt)]  = linprog(c_d(:,ii)+H_eq_d(:,:,ii)'*VV(:,ii),[],[],...
            G,g,LB_d(:,ii),UB_d(:,ii),options);
    
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

%% Animation
for i=1:N
    Task_assigned(i) =  find(ZZ_RA(:,i,end)>Threshold);
end

xx_actual_position = animation_EDE(xx_actual_position,Task_assigned,task,10);
fprintf('xx_actual_position . %g\n',xx_actual_position);
  
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
