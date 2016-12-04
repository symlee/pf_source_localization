function [P] = initializePF(nParticles, initDistribution, environmentSize)

P = zeros(nParticles, 2);

% TODO: might want to later spread out last line of points, dependent on
% shape of rectangle (taller or wider - currently fills up Y first)
nParticlesPerLine = ceil(sqrt(nParticles));

if isequal(initDistribution, 'uniform')
    dx = linspace(0, environmentSize(1), nParticlesPerLine + 2);
    dy = linspace(0, environmentSize(2), nParticlesPerLine + 2);
    [X, Y] = meshgrid(dx(2:length(dx)-1), dy(2:length(dy)-1));
    ind_x = reshape(X, [size(X,1)*size(X,2) 1]);
    ind_y = reshape(Y, [size(Y,1)*size(Y,2) 1]);

    P(1:nParticles,:) = [ind_x(1:nParticles) ind_y(1:nParticles)];
end

end