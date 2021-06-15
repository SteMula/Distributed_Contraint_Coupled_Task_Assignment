%
% DCS MATLAB Project 3, Task 2, Solution 3bis (Coupled with perturbations)
% Distributed Task Assignment for Robotic Networks
% Stefano Mulargia, Stefano Bortolotti
% April 2021
%
%%
% Starting from the coupled version we extended the problem also if there
% are perturbations.
%
%%
clear all; 
clc; 
close all;
addpath .\Functions;
rng(1)

% number of mobile robots: 2D (ground) agents and 3D (drones) agents
number_agent_ground = 5;
number_agent_drones = 5;

% N is the number of agents and number of tasks (>6 to avoid problems)
N = number_agent_ground + number_agent_drones;

% define the number of ground tasks. The remaining ones are in air
number_task_ground = 5; %randi(N,1,1);
number_task_drones = N - number_task_ground;
% If number_task_ground < number_agent_ground
% we will have more mobile robots in the same
% task (critical case)

% Number of iterations
MAXITERS = 1e4;

Threshold = 0.9;

% Note that if you decreasse the number of iteration you have to decrease
% also the threshold. Unfortunately there is not a formula to compute it
% but you have to make a trade off between number of iteration and the
% threshold.

% Options for linprog
options = optimoptions('linprog','Display','none');

% Initialization of task positions
task_2D = 10*rand(2,1,number_task_ground);
task_3D = 10*rand(3,1,number_task_drones);

% Initialization of agent positions
pos_init_mob_rob = 10*rand(2,1,number_agent_ground);
pos_init_drones = 10*rand(3,1,number_agent_drones);

%% Distance computation

% adding a third row of zeros
task_2D_copy = task_2D;
task_2D_copy(end+1,:,:) = 0;

% stack 2D tasks with 3D tasks coordinate matrices
task = zeros(3,1,N);
task(:,:,1:number_task_ground) = task_2D_copy;
task(:,:,number_task_ground+1:end) = task_3D;

% same procedure but with agents' coordinates
pos_init_mob_rob_copy = pos_init_mob_rob;
pos_init_mob_rob_copy(end+1, :,: ) = 0;
pos_init = zeros(3,1,N);
pos_init(:,:,1:number_agent_ground) = pos_init_mob_rob_copy;
pos_init(:,:,number_agent_ground+1:end) = pos_init_drones;

% compute agents-tasks distances
for i=1:N
    for j=1:N
    distance(i,j) = norm(pos_init(:,i)-task(:,j));
    end
end

% Perturbation vector Delta
Delta = rand(N,N);

% problem with perturbed cost
distance  = distance + Delta ;

% probability that a number in binord function will be 1. Namely we decide
% the total amount of task that all the agents are able to do.
p_S = 0.6;

% definition of feasible combinations: tasks (rows) and agents (columns)
% if 1 the agent can reach the task, if 0 it can't
% called S as "selection matrix"
S_mb = binornd(1,p_S,number_task_ground,number_agent_ground);
S_dr = binornd(1,p_S,N,number_agent_drones);

% adding rows of zeros in order to have rows=N 
S_zero = zeros(N-number_task_ground,number_agent_ground);
S_mb(end+1:end+N-number_task_ground,:) = S_zero;

% join mobile robots (2D) with drones (3D)
S = horzcat(S_mb,S_dr);

% only distances of feasible paths, the others are zeros
distance_mask = distance.*S;

% Since each z_i can be of different size we cannot use a matrix, so we
% decided to use a data structure called "cell"
c_distance = num2cell(distance_mask,1);

% From the cells we remove zeros
for i = 1:length(c_distance)
    c_distance{i}(c_distance{i}==0) = [];
end

%% Network graph
% p = 0.1;
% [AA_NW, AA] = binomialGraph(p, N, 'doubly')
[AA_NW, AA] = cyclicGraph(N, 'doubly');

%% Distributed Solution

% Once we computed the distance we can define the cost matrix "c_d", as
% the concatenation of the mobile robot distances from their tasks and the
% drones ones
c_d = {c_distance};
c_d = horzcat(c_d{:});

% Variables' initialization 
ZZ=cell(1,length(c_d),MAXITERS);

for i=1:MAXITERS
    for j=1:length(c_d)
        ZZ{1,j,i} = zeros(length(c_d{:,j}),1) ; 
    end
end
ZZ_RA = ZZ;

%%% Lambda
Lambda = zeros(N,N,MAXITERS); 

%%% Bounds
LB_d = ZZ(:,:,1);
UB_d = ZZ(:,:,1);
for j = 1:length(UB_d)
    UB_d{1,j} = ones(length(UB_d{:,j}),1); 
end

%%% Equality constraints

b_eq = ones(N)/N;
for ii = 1:N
    H_eq(:,:,ii) = diag(S(:,ii));
    [r_Heq(ii),~] = size(cell2mat(ZZ(:,ii,1)));
end
H_eq(:,all(~any(H_eq),1)) = [];

% H_eq definition
shift = 1;

for i=1:N
    H_eq_d{:,i} = H_eq(:,shift:shift+r_Heq(i)-1);
    shift = shift+r_Heq(i);
end

G_eq_d = H_eq_d;

for ii = 1:N
   G_eq_d{:,ii} = [repmat(0,N,length(ZZ{:,ii,1}))];
   G_eq_d{:,ii}(ii,:) = [repmat(1,1,length(ZZ{:,ii,1}))];
end

d = eye(N);
VV = zeros(N,N);
primal_cost = cell(MAXITERS,1);
dual_cost = cell(MAXITERS,1);
primal_cost_RA = cell(MAXITERS,1);
dual_cost_RA = cell(MAXITERS,1);
consensus_err = cell(MAXITERS,1);

for j = 1:length(primal_cost)
    primal_cost{j,1} = 0 ; 
    dual_cost{j,1} = 0 ; 
    primal_cost_RA{j,1} = 0 ; 
    dual_cost_RA{j,1} = 0 ; 
    consensus_err{j,1} = 0 ; 
end

gamma_const = 3;
gamma_exp = 0.6;

%% Centralized solution
c_d_length=0;
for ii=1:N
    c_d_length = c_d_length + length(c_d{1,ii});
end
c_centr = zeros(1,c_d_length);
jj = 1;
for ii=1:N
    c_centr(1,jj:jj+length(c_d{1,ii})-1) = cell2mat(c_d(1,ii))';
    jj = jj + length(c_d{1,ii});
end
G_centr =cell2mat(G_eq_d);
g_centr = ones(N,1);
H_eq_centr = H_eq;
b_eq_centr = ones(N,1);
UB_centr = ones(c_d_length,1);
LB_centr = zeros(c_d_length,1);

Equal = vertcat(H_eq_centr,G_centr);
b_equal = vertcat(b_eq_centr,g_centr);

[xopt, fopt] = linprog(c_centr,[],[],Equal,b_equal,LB_centr,UB_centr,options)

%% Distributed dual subgradient
for tt = 1:MAXITERS-1
  if mod(tt,100)==0
      fprintf('Iteration n. %d\n',tt);
  end

  gamma_t = gamma_const*(1/tt)^gamma_exp; % 1/tt^alpha with alpha in (0.5, 1]

  for ii=1:N
    N_ii = find(AA_NW(:,ii) == 1)';

    VV(:,ii) = AA(ii,ii) * Lambda(:,ii,tt);
    
    for kk = N_ii
      VV(:,ii) = VV(:,ii) + AA(ii,kk) * Lambda(:,kk,tt);
    end
  end

  % Primal Update

  for ii=1:N
    [ZZ{:,ii,tt}] = linprog(c_d{:,ii}+H_eq_d{:,ii}'*VV(:,ii),...
        [],[],G_eq_d{:,ii},...
        d(:,ii),LB_d{:,ii},UB_d{:,ii},options);
  end

  % Running average
  for ii=1:N
      if tt==1
          ZZ_RA(:,ii,tt) = ZZ(:,ii,tt);
      else
          ZZ_RA(:,ii,tt) = mat2cell((1/tt)*((tt-1)*cell2mat(ZZ_RA(:,ii,tt-1))+...
              cell2mat(ZZ(:,ii,tt))),length(ZZ_RA{:,ii,tt}),1);
      end
  end
  
  % Dual Update
  for ii=1:N
    grad_ii = H_eq_d{:,ii}*ZZ{:,ii,tt}- b_eq(:,ii);

    Lambda(:,ii,tt+1) = VV(:,ii) + gamma_t*grad_ii;
  end

  % Performance check
    Lambda_avg = mean(Lambda(:,:,tt),2);

  
  for ii=1:N
    c_d_length_ii = cell2mat(c_d(:,ii))' * cell2mat(ZZ(:,ii,tt));
    primal_cost(tt) = mat2cell(cell2mat(primal_cost(tt)) + c_d_length_ii,1,1);
    
    c_d_length_ii = cell2mat(c_d(:,ii))' * cell2mat(ZZ_RA(:,ii,tt));
    primal_cost_RA(tt) = mat2cell(cell2mat(primal_cost_RA(tt)) + c_d_length_ii,1,1);
    
    qq_ii = cell2mat(c_d(:,ii))' * cell2mat(ZZ(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ(:,ii,tt))-b_eq(:,ii));
    
    dual_cost(tt) =  mat2cell(cell2mat(dual_cost(tt)) + qq_ii,1,1);
    
    qq_ii = cell2mat(c_d(:,ii))' * cell2mat(ZZ_RA(:,ii,tt)) + ...
        Lambda(:,ii,tt)' * (cell2mat(H_eq_d(:,ii)) * cell2mat(ZZ_RA(:,ii,tt))-b_eq(:,ii));

    dual_cost_RA(tt) = mat2cell(cell2mat(dual_cost_RA(tt)) + qq_ii,1,1);

    consensus_err(tt) = mat2cell(cell2mat(consensus_err(tt)) + ...
        norm(Lambda(:,ii,tt) - Lambda_avg),1,1);
  end
end
%% Last value [for plot]

tt = MAXITERS;
fprintf('Iteration n. %d\n',tt);
	Lambda_avg = mean(Lambda(:,:,tt),2);

for ii=1:N
    N_ii = find(AA_NW(:,ii) == 1)';

    VV(:,ii) = AA(ii,ii)*Lambda(:,ii,tt);

    for kk = N_ii
        VV(:,ii) = VV(:,ii) + AA(ii,kk)*Lambda(:,kk,tt);
    end
end

for ii=1:N

    for jj=1:N
        [ZZ{:,jj,tt}] = linprog(c_d{:,jj}+H_eq_d{:,jj}'*VV(:,jj),...
            [],[],G_eq_d{:,jj},...
            d(:,jj),LB_d{:,jj},UB_d{:,jj},options);
    end

    ZZ_RA(:,ii,tt) = mat2cell((1/tt)*((tt-1)*cell2mat(ZZ_RA(:,ii,tt-1))+...
        cell2mat(ZZ(:,ii,tt))),length(ZZ_RA{:,ii,tt}),1);

    c_d_length_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ(:,ii,tt));
    primal_cost(tt) = mat2cell(cell2mat(primal_cost(tt)) + c_d_length_ii,1,1);

    c_d_length_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ_RA(:,ii,tt));
    primal_cost_RA(tt) = mat2cell(cell2mat(primal_cost_RA(tt)) + c_d_length_ii,1,1);

    qq_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ(:,ii,tt))-b_eq(:,ii));

    dual_cost(tt) =  mat2cell(cell2mat(dual_cost(tt)) + qq_ii,1,1);

    qq_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ_RA(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ_RA(:,ii,tt))-b_eq(:,ii));

    dual_cost_RA(tt) = mat2cell(cell2mat(dual_cost_RA(tt)) + qq_ii,1,1);

end

%% task assignment

% find the right one inside every column of z (so the feasible task)
for i = 1:N
    Task_assigned_index(i) = find(cell2mat(ZZ_RA(:,i,end))>Threshold);
end

% assign the tasks to the agents (the elements of Task_augment are tasks,
% the columns are the agents)
for ii=1:N
   [r,c] = find(S(:,ii));
   Task_augment(ii) = r(Task_assigned_index(ii));
end

%% Animation 
animation_task2(pos_init,Task_augment,task,N);
fprintf('actual position . %d\n',pos_init);

%% plot dual cost
figure
  semilogy(1:MAXITERS,abs(cell2mat(dual_cost(1:MAXITERS))-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogy(1:MAXITERS,abs(cell2mat(dual_cost_RA(1:MAXITERS))-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('dual cost error')
  legend('dual cost','dual cost with running avg')

%% plot primal cost
figure
  semilogx(1:MAXITERS,abs(cell2mat(primal_cost(1:MAXITERS))-fopt), 'LineWidth',2);
  hold on, grid on, zoom on
  semilogx(1:MAXITERS,abs(cell2mat(primal_cost_RA(1:MAXITERS))-fopt), 'LineWidth',2);
  xlabel('t')
  ylabel('primal cost error')
  legend('primal cost','primal cost with running avg')

%% plot primal vs dual cost
% figure
%   semilogy(1:MAXITERS,abs(primal_cost(1:MAXITERS)-fopt), 'LineWidth',2);
%   hold on, grid on, zoom on
%   semilogy(1:MAXITERS,abs(primal_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
%   semilogy(1:MAXITERS,abs(dual_cost(1:MAXITERS)-fopt), 'LineWidth',2);
%   semilogy(1:MAXITERS,abs(dual_cost_RA(1:MAXITERS)-fopt), 'LineWidth',2);
%   xlabel('t')
%   ylabel('cost error')
%   legend('primal cost','primal cost with running avg','dual cost','dual cost with running avg')
%   
%% plot lambda
Lambdax=reshape(Lambda,N*N,tt);

figure
  plot(1:MAXITERS,Lambdax(:,:), 'LineWidth',1.5);
  hold on, grid on, zoom on
  xlabel('t')
  ylabel('\lambda_i^t')

%% plot consensus error
figure
  semilogy(1:MAXITERS,cell2mat(consensus_err(1:MAXITERS)), 'LineWidth',2);
  hold on, grid on, zoom on
  refresh
  xlabel('t')
  ylabel('consensus error')
  
%%
% mean(abs(cell2mat(primal_cost(1:MAXITERS))-fopt))
% mean(abs(cell2mat(dual_cost(1:MAXITERS-2))-fopt))
% mean(cell2mat(consensus_err(:)))

%% comparison between nominal and perturbated versions
% first run the coupled version and import workspace. Keep only consesus_err and
% dual_cost and rename them as consensus_err_norm and dual_cost_norm. Then
% run this version and rename consesus_err and dual_cost as consensus_err_pert and
% dual_cost_pert. At the end, plot.

%%
% clearvars -except  fopt consensus_err_pert dual_cost_pert dual_cost_norm consensus_err_norm MAXITERS fopt_norm
% 
% %% comparison consensus error
% figure
%   semilogy(1:MAXITERS,cell2mat(consensus_err_norm(1:MAXITERS)), 'LineWidth',2);
%   hold on, grid on, zoom on
%   semilogy(1:MAXITERS,cell2mat(consensus_err_pert(1:MAXITERS)), 'LineWidth',2);
%   refresh
%   legend('consensus error','consensus error perturbated')
%   xlabel('t')
%   ylabel('consensus error')
%   
% %% comparison dual cost
% figure
%   semilogy(1:MAXITERS-2,abs(cell2mat(dual_cost_norm(1:MAXITERS-2))-fopt_norm), 'LineWidth',2);
%   hold on, grid on, zoom on
%   semilogy(1:MAXITERS-2,abs(cell2mat(dual_cost_pert(1:MAXITERS-2))-fopt_pert), 'LineWidth',2);
%   xlabel('t')
%   ylabel('cost error')
%   legend('dual cost','dual cost perturbated')
%   