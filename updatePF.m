function [ pf_t ] = updatePF(pf_t_1, z_t, sources, Xr)
%UPDATEPF Summary of this function goes here
%   Detailed explanation goes here

nParticles = size(pf_t_1, 1);

[X_t_hat, X_t] = deal(zeros(size(pf_t_1)));
w_t = zeros(size(pf_t_1, 1), 1);

% sampling
for particle = 1:nParticles
    x_t_1 = pf_t_1(particle,:);
    x_t = x_t_1;
    X_t_hat(particle, :) = x_t;
    w_t(particle) = weight(z_t, x_t, sources, Xr);
end

% normalize weights and find cdf
w_t = w_t ./ sum(w_t);
w_t_cum = cumsum(w_t);
% sum(w_t)

% calculate effective number of particles (estimates efficiency of
% representation)
Neff = 1 / sum(w_t.^2);

% resampling according to weights
X_t = resample(X_t_hat, w_t_cum, X_t);

pf_t = X_t;

end

function [X_t] = resample(X_t_hat, w_t_cum, X_t)
    for particle = 1:size(X_t,1)
        X_t(particle,:) = X_t_hat(find(rand <= w_t_cum, 1),:);
    end
end

% affected - need to sample divergence from all known number of sources
function [w_t] = weight(z_t, x_t, sources, Xr)
    measurementExpected = 0;
    
    for j = 1:length(sources)
        measurementExpected = measurementExpected + sources{j}.str / (10 ^ (sources{j}.n * log10(norm(Xr - x_t))));
    end
    
%     % perfect sensor model
%     if abs(measurementExpected - z_t) < 0.001
%         w_t = 1;
%     else 
%         w_t = 0;
%     end

    % need to normalize divergence to work with weak measurements
    divergence = (measurementExpected - z_t)/ z_t;

    % highest weight given to zero divergence between noisy measurement and
    % expected measurement given expected source location and noisy position 
    w_t = normpdf(divergence, 0, 0.2);
end