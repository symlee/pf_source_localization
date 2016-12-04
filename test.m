% dt = 0.05; % time step
% initialPose = [0  0  0  0]';
% carbot = ExampleHelperCarBot(initialPose, dt); 
% 
% simulationTime = 0;
% for simulationtime = 1:20
%     uCmd(1) = 0.7*abs(sin(simulationTime)) + 0.1;  % linear velocity
%     uCmd(2) = 0.08*cos(simulationTime);            % angular velocity
% 
%     drive(carbot, uCmd);
%     updatePlot(carbot)
% end

%% drawing from probability distribution

A = [1 0 5]';
B = A ./ sum(A);
C = cumsum(B)

count = zeros(3,1);

for j = 1:1000000
    ind = find(rand <= C,1);
    count(ind,1) = count(ind) + 1;
end

count

%% Looking at Normal Distributions of Mobility and Sensor models

x = [-5:0.1:5];
norm = normpdf(x, 0, 0.5);
figure(1)
plot(x, norm)

