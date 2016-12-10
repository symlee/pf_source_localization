close all
clear all

% parameters of environment
nGridLengthX = 100;         % determines number of sides along grid
nGridLengthY = 100;         % determines number of sides along grid
grid = zeros(nGridLengthY, nGridLengthX);
circleSize = 108;

% parameters of source
nSources = 1;
sources = cell(nSources, 1);
for j = 1:nSources
    source.str  = ((nGridLengthX + nGridLengthY)/2) ^ 2;                    % source strength
    source.n    = 3;                      % decay exponent (for free space)
    source.loc  = [nGridLengthX/2, nGridLengthY/2];               % location of source (x, y)

    sources{j} = source;
end

% parameters for robot
Xr = [1, 1];                % location of robot
stepSize = nGridLengthX * 0.2;
varMobility = 0.00002;
varSensor = 0.00002;

nParticles = 1000;
initDistribution = 'uniform';

% initialize particle filter
pf_t_1 = initializePF(nSources, nParticles, initDistribution, [nGridLengthX, nGridLengthY]);

figure(1)
hold on
scatter(pf_t_1(:,1), pf_t_1(:,2))
scatter(Xr(1), Xr(2), circleSize, 'ro', 'filled')
for j = 1:nSources
    scatter(sources{j}.loc(1), sources{j}.loc(2), circleSize, 'go', 'filled')
end
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
        
        % check validity of generated trajectory
        Xr_proposed = Xr + u_t;
        if (Xr_proposed(1) >= 0 && Xr_proposed(1) <= nGridLengthX && Xr_proposed(2) >= 0 && Xr_proposed(2) <= nGridLengthY)
            trajValid = 1;
              Xr_gt = Xr + u_t;

%               % perfect motion model
%               Xr = Xr_gt;
            
            % noisy motion model - sample from gaussian 
            Xr = normrnd(Xr_gt, varMobility);
        end
    end
    
    % receive measurement from where robot actually is but corrupted from
    % sensor noise
    z_t = 0;
    for j = 1:nSources
        nonCorruptedReading = sources{j}.str / (10 ^ (sources{j}.n * log10(norm(sources{j}.loc - Xr_gt))));
        corruptedReading = abs(normrnd(nonCorruptedReading, varSensor));
        z_t = z_t + corruptedReading;
    end
    
    pf_t = updatePF(pf_t_1, z_t, sources, Xr);

    % plot robot and particles
    figure(1)
    clf
    hold on
    scatter(pf_t(:,1), pf_t(:,2))
    scatter(Xr(1), Xr(2), circleSize, 'ro', 'filled')
    for j = 1:nSources
        scatter(sources{j}.loc(1), sources{j}.loc(2), circleSize, 'go', 'filled')
    end

    axis([0 nGridLengthX 0 nGridLengthY])
    axis equal
    
    pf_t_1 = pf_t;
    iterations = iterations + 1;
end