function animation_task2(xx_actual_position,Task_assigned,task,numberofpartition)
% Basically you pass to this function the actual posizion of the agent,
% then you tell him the destination. With linpace you interpolate the two
% point, difine the trajectory as a simple line, then plot it


%   xx_actual_position  is a matrix NxN collecting every position of the
%                       agent in the space 
%   Task_assigned       tell you which task is assigne i.e A1 has Task 2,
%                       so we have 2  
%   task                is a matrix NxN collecting every position of the
%                       task in the space 
%   numberofpartition   is the number of partition that linspace will use


[N,M] = size(xx_actual_position) ;                                             % xx_actual_position is squared
                                                                           % so the number of col and row are the same
xx_actual_position = reshape(xx_actual_position,3,M);
task = reshape(task,3,M);
if numberofpartition > 5
    numberofpartition= 10;
end    

nPoints = numberofpartition ;%number of points in each row

x1 = linspace(0,1,nPoints);     % define a set of point 


% A is NxN matrix initialized to zeros collecting the future point
% obtainter by creating a line in the space starting from the actual
% postion of each agent and the task position the they have to do

A = zeros(N,nPoints,M); % position, number of points, number of agent
for i=1:M
    A(:,:,i) = xx_actual_position(:,i) + x1.*(task(:,Task_assigned(i))- xx_actual_position(:,i));
end
for i=1:M
    B(:,:,i) =A(:,:,i)'; % number of agent, number of points, position
end
%% Animation
grid on;
for i=1:M   %animation robot
h(i)= animatedline('Color','r','MaximumNumPoints', 10,'Marker','o');
end
for i=1:M   %animation task
j(i)= animatedline('Color','b','MaximumNumPoints', 10,'Marker','o');
end

% I'm assuming that all coordinates are 3-dimensional, you can just set to
% zero the third one 
stopiter = 10;
x = B(:,1,:);
y = B(:,2,:);
if N == 3
    z = B(:,3,:);
else
    z = zeros(stopiter,nPoints,1);
end


% Force a 3D view
view(3);
% I set the following for-cycle to 2 because I want that the agent do only
% one step at each iteration
% for i=1:N
for k = 1:stopiter
    addpoints(j(k),task(1,k),task(2,k),task(3,k));
    drawnow
end
for k = 1:stopiter
    for i=1:stopiter
    addpoints(h(i),x(k,:,i),y(k,:,i),z(k,:,i));     
    end
    drawnow
    pause(1)
end

end

