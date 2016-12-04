close all
clear all

nGridLengthX = 100;         % determines number of sides along grid
nGridLengthY = 100;         % determines number of sides along grid
grid = zeros(nGridLengthY, nGridLengthX);

P0 = 100;                   % source strength
n = 2;                      % decay exponent (for free space)
Xs = [50, 50];              % location of source (x, y)
Xr = [1, 1];                % location of robot
stepSize = nGridLengthX * 0.1;

save('field.mat', 'P0', 'n', 'Xs')

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
    goDown = randi([0, 1]);
    goRight = randi([0, 1]);

    % randomly explore map
    if goDown == 1
        if goRight == 1
            u_t = [stepSize -stepSize];     % (x, y)
        else
            u_t = [-stepSize -stepSize];
        end
    else
        if goRight == 1
            u_t = [stepSize stepSize];
        else
            u_t = [-stepSize stepSize];
        end
    end
    Xr = Xr + u_t;
    
    % receive measurement
    z_t = P0 - 10 * n * log10(norm(Xs - Xr));
    
    pf_t = updatePF(pf_t_1, u_t, z_t, Xr);
    
    % plot robot and particles
    figure(1)
    clf
    hold on
    scatter(pf_t(:,1), pf_t(:,2))
    scatter(Xr(1), Xr(2), 'r+')
    scatter(Xs(1), Xs(2), 'go')
    
    pf_t_1 = pf_t;
    iterations = iterations + 1;
end