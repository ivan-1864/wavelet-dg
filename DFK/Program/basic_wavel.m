function [psi] = basic_wavel(t)
% create basic wavelet

% mhat wavelet
psi = (1-t.^2).*exp(-t.^2/2);

% % haar wavelet
% len = length(t);
% psi = zeros(len, 1) ...
%       + ones(len, 1) .* (t > 0) .* (t < 1/2) ...
%       - ones(len, 1) .* (t >= 1/2) .* (t < 1);
end