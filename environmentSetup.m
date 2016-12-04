close all
clear all

global Xs

% parameters of environment
nGridLengthX = 100;         % determines number of sides along grid
nGridLengthY = 100;         % determines number of sides along grid
grid = zeros(nGridLengthY, nGridLengthX);

% parameters of source
P0 = 100;                   % source strength
n = 2;                      % decay exponent (for free space)
Xs = [50, 50];              % location of source (x, y)

% parameters for robot
Xr = [1, 1];                % location of robot
stepSize = nGridLengthX * 0.1;
varMobility = 0.2;
varSensor = 0.2;

% save('field.mat', 'P0', 'n', 'Xs');
% save('mobility_sensor.mat', 'varMobility', 'varSensor');
nParticles = 4900;
initDistribution = 'uniform';

% initialize particle filter
pf_t_1 = initializePF(nParticles, initDistribution, [nGridLengthX, nGridLengthY]);

figure(1)
hold on
scatter(pf_t_1(:,1), pf_t_1(:,2))
scatter(Xr(1), Xr(2), 'r+')
scatter(Xs(1), Xs(2), 'go')
axis([0 nGridLengthX 0 nGridLengthY])
axis equal

converged = 0;
iterations = 0;
while converged == 0
    % randomly explore map
    trajValid = 0;
    while ~trajValid
        X_proposed = [0 0];
        motion = randi([0, 3]);
        
        switch motion
            case 0
                u_t = [0 stepSize];
            case 1
                u_t = [0 -stepSize];
            case 2
                u_t = [stepSize 0];
            case 3
                u_t = [-stepSize 0];
        end

        Xr_proposed = Xr + u_t;
        if (Xr_proposed(1) >= 0 && Xr_proposed(1) <= nGridLengthX && Xr_proposed(2) >= 0 && Xr_proposed(2) <= nGridLengthY)
            trajValid = 1;
            Xr = Xr + u_t;
        end
    end
    
    % receive measurement
    z_t = P0 - 10 * n * log10(norm(Xs - Xr));
    
    pf_t = updatePF(pf_t_1, u_t, z_t, P0, n, Xr, varMobility, varSensor);
    
    % plot robot and particles
    figure(1)
    clf
    hold on
    scatter(pf_t(:,1), pf_t(:,2))
    scatter(Xr(1), Xr(2), 'r+')
    scatter(Xs(1), Xs(2), 'go')
    axis([0 nGridLengthX 0 nGridLengthY])
    axis equal
    
    pf_t_1 = pf_t;
    iterations = iterations + 1;
end