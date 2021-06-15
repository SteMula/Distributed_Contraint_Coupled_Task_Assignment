%
% DCS MATLAB Project 3, Task 2, Solution 2 (Decoupled)
% Distributed Task Assignment for Robotic Networks
% Stefano Mulargia, Stefano Bortolotti
% April 2021
%
%%
% In this code we developed a decoupled task assignment in which ground
% tasks can only be assigned to mobile robots and air tasks can only be
% assigned to drones.
%
%%
clear all; 
clc; 
close all;
addpath .\Functions;
rng(1)

% N is the number of Agents and number of tasks (>6 to avoid problems)
N = 10;          

% define the number of task mobile robots and the remain ones are left for
% drones
number_task_ground = 5;
number_task_drones = N - number_task_ground;

% For simplicity we consider always that the number of tasks on the ground
% are equal to the number of mobile robots. The same holds for drones.
number_agent_ground = number_task_ground;
number_agent_drones = number_task_drones; 

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
pos_init_mob_rob = 10*rand(2,1,number_task_ground);
pos_init_drones = 10*rand(3,1,number_task_drones);

%% Distance computation

% compute the distance for mobile robots
for i=1:number_task_ground
    for j=1:number_task_ground
        distance_mb(i,j) = norm(pos_init_mob_rob(:,i)-task_2D(:,j));
    end
end

% Since some agent, for any reason, may be not allowed to compute some task
% we generate a random matrix S_mb that tell us which task an agent is
% allowed to do
p_S = 0.6;
S_mb = binornd(1,p_S,number_task_ground,number_task_ground);
distance_mb_mask = distance_mb.*S_mb;

% drones
for i=1:number_task_drones
    for j=1:number_task_drones
        distance_dr(i,j) = norm(pos_init_drones(:,i)-task_3D(:,j));
    end
end
S_dr = binornd(1,p_S,number_task_drones,number_task_drones);
distance_dr_mask = distance_dr.*S_dr;

% Since each z_i can be of different size we cannot use a matrix, the
% choice we made is to use a data structure called "cell"
c_distance_dr = num2cell(distance_dr_mask,1);

% From the cells we remove zeros (drones)
for i = 1:length(c_distance_dr)
    c_distance_dr{i}(c_distance_dr{i}==0) = [];
end

c_distance_mb = num2cell(distance_mb_mask,1);

% From the cells we remove zeros (mobile robots)
for i = 1:length(c_distance_mb)
    c_distance_mb{i}(c_distance_mb{i}==0) = [];
end

%% Network graph
p = 0.2;
[AA_NW, AA] = binomialGraph(p, N, 'doubly');
% [AA_NW, AA] = cyclicGraph(N, 'doubly');

%% Distributed Solution

% Once we have compute the distance we cen define the cost matrix "c_d", as
% the concatenation of the mobile robot distance from their task and the
% drones ones
c_d = {c_distance_mb,c_distance_dr};
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
for j=1:length(UB_d)
    UB_d{1,j} = ones(length(UB_d{:,j}),1) ; 
end

%%% Equality constraints
b_eq = ones(N)/N;
S = blkdiag(S_mb,S_dr);
for ii=1:N
    H_eq(:,:,ii)=diag(S     (:,ii));
end
H_eq( :, all( ~any( H_eq ), 1 ) ) = [];

for ii=1:N
 [rrr(ii),ccc(ii)] = size(cell2mat(ZZ(:,ii,1)));
end

% H_eq definition
jjj=1;

for i=1:N
    H_eq_d{:,i} = H_eq(:,jjj:jjj+rrr(i)-1);
    jjj = jjj+rrr(i);
end

G_eq_d = H_eq_d;

for ii=1:N
   G_eq_d{:,ii}=[repmat(0,N,length(ZZ{:,ii,1}))];
   G_eq_d{:,ii}(ii,:)=[repmat(1,1,length(ZZ{:,ii,1}))];
end

d=eye(N);

VV = zeros(N,N);
primal_cost = cell(MAXITERS,1);
dual_cost = cell(MAXITERS,1);
primal_cost_RA = cell(MAXITERS,1);
dual_cost_RA = cell(MAXITERS,1);
consensus_err = cell(MAXITERS,1);

for j=1:length(primal_cost)
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

% reshape for linprog
G_centr =cell2mat(G_eq_d);
g_centr = ones(N,1);
H_eq_centr = H_eq;
b_eq_centr = ones(N,1);
UB_centr = ones(c_d_length,1);
LB_centr = zeros(c_d_length,1);

Equal = vertcat(H_eq_centr,G_centr);
b_equal =  vertcat(b_eq_centr,g_centr);

[xopt, fopt] = linprog(c_centr,[],[],Equal,b_equal,LB_centr,UB_centr)

%% Distributed dual subgradient
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
        [ZZ{:,ii,tt}] = linprog(c_d{:,ii}+H_eq_d{:,ii}'*VV(:,ii),...
            [],[],nonzeros(cell2mat(G_eq_d(:,ii)))',...
            1,LB_d{:,ii},UB_d{:,ii},options);
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

  % Computation of primal cost, dual cost and consensus error
  for ii=1:N
    ff_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ(:,ii,tt));
    primal_cost(tt) = mat2cell(cell2mat(primal_cost(tt)) + ff_ii,1,1);
    
    ff_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ_RA(:,ii,tt));
    primal_cost_RA(tt) = mat2cell(cell2mat(primal_cost_RA(tt)) + ff_ii,1,1);
    
    qq_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ(:,ii,tt))-b_eq(:,ii));
    
    dual_cost(tt) =  mat2cell(cell2mat(dual_cost(tt)) + qq_ii,1,1);
    
    qq_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ_RA(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ_RA(:,ii,tt))-b_eq(:,ii));

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

    % v_i =  sum_{k \in N_i U { i } } w_ik mu_k^t
    VV(:,ii) = AA(ii,ii)*Lambda(:,ii,tt);
    
    for kk = N_ii
    VV(:,ii) = VV(:,ii) + AA(ii,kk)*Lambda(:,kk,tt);
    end
  end
  
for ii=1:N
    
  for jj=1:N 
        [ZZ{:,jj,tt}] = linprog(c_d{:,jj}+H_eq_d{:,jj}'*VV(:,jj),...
            [],[],nonzeros(cell2mat(G_eq_d(:,jj)))',...
            1,LB_d{:,jj},UB_d{:,jj},options);
  end
                
    ZZ_RA(:,ii,tt) = mat2cell((1/tt)*((tt-1)*cell2mat(ZZ_RA(:,ii,tt-1))+...
              cell2mat(ZZ(:,ii,tt))),length(ZZ_RA{:,ii,tt}),1);
          
    ff_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ(:,ii,tt));
    primal_cost(tt) = mat2cell(cell2mat(primal_cost(tt)) + ff_ii,1,1);
    
    ff_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ_RA(:,ii,tt));
    primal_cost_RA(tt) = mat2cell(cell2mat(primal_cost_RA(tt)) + ff_ii,1,1);
    
    qq_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ(:,ii,tt))-b_eq(:,ii));
    
    dual_cost(tt) =  mat2cell(cell2mat(dual_cost(tt)) + qq_ii,1,1);
    
    qq_ii = cell2mat(c_d(:,ii))'*cell2mat(ZZ_RA(:,ii,tt)) + ...
        Lambda(:,ii,tt)'*(cell2mat(H_eq_d(:,ii))*cell2mat(ZZ_RA(:,ii,tt))-b_eq(:,ii));

    dual_cost_RA(tt) = mat2cell(cell2mat(dual_cost_RA(tt)) + qq_ii,1,1);

end

%%
pos_init_mob_rob_copy=pos_init_mob_rob;
pos_init_mob_rob_copy( end+1, :,: ) = 0; 
task_2D_copy=task_2D;
task_2D_copy( end+1, :,: ) = 0; 
xx_actual_position =zeros(3,1,10);
xx_actual_position(:,:,1:number_agent_ground)=pos_init_mob_rob_copy;
xx_actual_position(:,:,number_agent_ground+1:end)=pos_init_drones;

task =zeros(3,1,10);
task(:,:,1:number_agent_ground)=task_2D_copy;
task(:,:,number_agent_ground+1:end)=task_3D;
for i=1:N
    Task_assigned_index(i) =  find(cell2mat(ZZ_RA(:,i,end))>Threshold);
end

%% task assignment
for ii=1:number_agent_ground
   [r,c] = find(S_mb(:,ii));
   Task_augment(ii) = r(Task_assigned_index(ii));
end

for ii=1:number_agent_drones
   [r,c] = find(S_dr(:,ii));
   Task_augment(ii+number_agent_ground) = r(Task_assigned_index(ii+number_agent_ground))+number_agent_ground;
end

%% Animation 

animation_task2(xx_actual_position,Task_augment,task,10);

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
  semilogy(1:MAXITERS,abs(cell2mat(primal_cost(1:MAXITERS))-fopt), 'LineWidth',2);
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