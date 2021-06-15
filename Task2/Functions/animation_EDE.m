function [xx_actual_position] = animation_EDE(xx_actual_position,Task_assigned,task,numberofpartition)
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


[N,M] = size(xx_actual_position)                                              % xx_actual_position is squared
                                                                           % so the number of col and row are the same

if numberofpartition > 5
    numberofpartition= 10;
end    

nPoints = numberofpartition %number of points in each row

x1 = linspace(0,1,nPoints);     % define a set of point 


% A is NxN matrix initialized to zeros collecting the future point
% obtainter by creating a line in the space starting from the actual
% postion of each agent and the task position the they have to do

A = zeros(N,nPoints,N); % position, number of points, number of agent
for i=1:N
    A(:,:,i) = xx_actual_position(:,i) + x1.*(task(:,Task_assigned(i))- xx_actual_position(:,i));
end
for i=1:N
    B(:,:,i) =A(:,:,i)'; % number of agent, number of points, position
end
%% Animation
grid on;
h1 = animatedline('Color','b','MaximumNumPoints', 10,'Marker','o');
h2 = animatedline('Color','r','MaximumNumPoints', 10,'Marker','o');
h3 = animatedline('Color','g','MaximumNumPoints', 10,'Marker','o');
h4 = animatedline('Color','k','MaximumNumPoints', 10,'Marker','o');
h5 = animatedline('Color','k','MaximumNumPoints', 10,'Marker','o');
h6 = animatedline('Color','k','MaximumNumPoints', 10,'Marker','o');

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

B
x

% Force a 3D view
view(3);
% I set the following for-cycle to 2 because I want that the agent do only
% one step at each iteration
% for i=1:N
    for k = 1:stopiter
        addpoints(h1,x(k,:,1),y(k,:,1),z(k,:,1)); %,task(:,1),task(:,2),task(:,3));        
        addpoints(h2,x(k,:,2),y(k,:,2),z(k,:,2));
        addpoints(h3,x(k,:,3),y(k,:,3),z(k,:,3));
%         addpoints(h4,task(1,1),task(2,1),task(3,1));
%         addpoints(h5,task(2,1),task(2,2),task(2,3));
%         addpoints(h6,task(3,1),task(3,2),task(3,3));
        addpoints(h4,task(1,1),task(2,1),task(3,1));
        addpoints(h5,task(1,2),task(2,2),task(3,2));
        addpoints(h6,task(1,3),task(2,3),task(3,3));
        drawnow
        pause(1)
    end
% end
xx_actual_position = [x(2,:,1),x(2,:,2),x(2,:,3);
                        y(2,:,1),y(2,:,2),y(2,:,3);
                        z(2,:,1),z(2,:,2),z(2,:,3)]
end

