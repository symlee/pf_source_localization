function [ pf_t ] = updatePF( pf_t_1, u_t, z_t, Xr )
%UPDATEPF Summary of this function goes here
%   Detailed explanation goes here

load('field.mat');

nParticles = size(pf_t_1, 1);

% [X_t, X_t_hat] = deal([]);
[X_t_hat, X_t] = deal(zeros(size(pf_t_1)));
w_t = zeros(size(pf_t_1, 1), 1);

% sampling
for particle = 1:nParticles
    x_t_1 = pf_t_1(particle,:);
    x_t = sample(u_t, x_t_1);
    X_t_hat(particle, :) = x_t;
    w_t(particle) = weight(z_t, x_t, P0, n, Xr);
end

% normalize weights and find cdf
w_t = w_t ./ sum(w_t);
w_t_cum = cumsum(w_t);
sum(w_t)
% resampling according to weights
for particle = 1:nParticles
    X_t(particle,:) = X_t_hat(find(rand <= w_t_cum, 1),:);
end

pf_t = X_t;

end
function [x_t] = sample(u_t, x_t_1)
    % mobility model is currently perfect, so can only draw from a
    % single sample with certainty 
    x_t = x_t_1 - u_t;  
    
%     maybe sample from gaussian centerd around x_t_1 - u_t
end

function [w_t] = weight(z_t, x_t, P0, n, Xr)
    %  sample from gaussian centered around the one ring where you can get a good measurement
    if abs(P0 - 10 * n * log10(norm(Xr - x_t)) - z_t) < 0.2
        w_t = 1;
    else 
        w_t = 0;
    end
end