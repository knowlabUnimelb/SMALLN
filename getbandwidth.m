function h = getbandwidth(data)
% Turner2014 - Eq. 14

N = numel(data);
% h = .9 * min(std(data), (prctile(data, 75) - prctile(data,25))./1.34) *
% N.^(-1/5); % Silvermans rule of thumb

I = iqr(data);
S = std(data);
h = .9 * min(S, I) * N.^(-1/5); % Silvermans rule of thumb
