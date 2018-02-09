% 2017 Small-N project - Additive factors simulation (Smith & Little)
% Simulate small number of subjects with probability of interaction set to p
% Analyse individuals using logNormal: use both weighted least squares and fit both full and constrained (i.e., no interaction) lognormal regression models
% Analyse groups using ANOVA
%
% This code generates a large nubmer of iterations of small subjects
% comparing individual vs group analyses
%
% Increase sample size to group levels in powers of 2

%%
clear all; clc; close all force hidden

%% Stimulation parameters
n = 128;    % Number of subjects (small n)
M = 500;  % Data points per subject
nBootstrap = 1000;
nbootsamples = 400;
N = 100; % Number of simulation iterations (on each iteration a new set of subjects is generated)

pset = [.1 .25 .5 .75 .9];

%% Fixed GLM parameters for each subject
factor1 = [0 1]; % Levels of factor 1
factor2 = [0 1]; % Levels of factor 2
con = allcomb(factor1, factor2); % All experimental conditions

intercept   = 400;       % Baseline RT
beta        = [150 200]; % Weight of factor 1 and factor 2 (main effects)
interaction_mean = 50;  % Mean of the interaction distribution
interaction_std = 5;    % Std of the interaction distribution (i.e., it is reasonable to assume that some proportion of subjects won't show the effec)
v           = 100^2;     % Variance on the RT scale - Increased variance results in more skew which shows how the least squares method on RT scale fails

%% Conversion functions for Lognormal model for MLE
% If X ~ Lognormal(mu, sigma), then log(X) ~ Normal(mu, sigma)
% m, s are the arithmetic mean and variance of X (on non-log scale) --> RT scale
% mu, sigma are parameters of the logNormal
%   In general, we want to apply the regression on the non-log scale to estimate the
%   parameters of main effects and interaction. This requires converting the additive model to mu and sigma
% Various conversions:
mufun    = @(m, v)(log((m.^2)./sqrt(v+m.^2)));
mufun2   = @(m, sigma)((log(m) - (sigma.^2)/2)); % e.g., Gustavsson2014, p. 3523-3524.
sigmafun = @(m, v)(sqrt(log(1+ (v./(m.^2)))));
mfun     = @(mu, v)(2^(1/2)*v^(1/2).*exp(mu./2).*(-1./(exp(mu) - (4*v + exp(2*mu)).^(1/2))).^(1/2)); % Convert mu and v to m (e.g., from log to non-log mean)
mfun2    = @(mu, sigma)(exp(sigma.^2/2 + mu));
vfun     = @(m, sigma)(m.^2.*(exp(sigma.^2) - 1));

%% Functions for simulating data
simmeans = @(intercept, beta)([intercept; intercept + beta(2); intercept + beta(1); intercept + beta(1) + beta(2) + beta(3)]);
simmeans_con = @(intercept, beta)([intercept; intercept + beta(2); intercept + beta(1)]);

%% Lognormal models for parameter estimation (negative log likelihoods)
% Full model: X includes intercept, factor 1, factor 2, and interaction
lnL     = @(parms, y,  X)(-sum(-log(y)...
    -.5 * log(2 * pi * sigmafun(X * parms(1:4), parms(5).^2).^2)...
    -.5 * (1./(sigmafun(X * parms(1:4), parms(5).^2)).^2) .* (log(y) - mufun(X * parms(1:4), parms(5).^2)).^2));

% Constrained model: X includes intercept, factor 1, factor 2 but no interaction
lnL_con = @(parms, y,  X)(-sum(-log(y)...
    - .5 * log(2 * pi * sigmafun(X * parms(1:3), parms(4).^2).^2)...
    - .5 * (1./(sigmafun(X * parms(1:3), parms(4).^2).^2)) .* (log(y) - mufun(X * parms(1:3), parms(4).^2)).^2));

%% Options for fminsearch
defopts = optimset ('fminsearch');
options = optimset (defopts, 'Display', 'off', 'MaxFunEvals', 1e5, 'MaxIter', 1e5, 'TolFun', 1e-6);

%% Output table
% Column headers for each analysis
indivResults(1,:) = {'Subject', 'Intercept', 'Beta_1', 'Beta_2', 'Beta_1x2', 'std',...
    'wls_int',   'wls_b1',  'wls_b2', 'wls_b1x2', 'wls_p_int', 'wls_p_b1',  'wls_p_b2', 'wls_p_b1x2', 'sigma_est_wls',...
    'full_int',  'full_b1', 'full_b2', 'full_b1x2', 'full_sigma', 'con_int',  'con_b1', 'con_b2', 'con_sigma',...
    'G2', 'pval'};
groupResults(1,:) = {'S1_int', 'S2_int', 'S3_int', 'S4_int', 'df_btw_A', 'df_wth_A', 'F_A', 'p_A',...
    'df_btw_B', 'df_wth_B', 'F_B', 'p_B', 'df_btw_A*B', 'df_wth_A*B', 'F_A*B', 'p_A*B'};

subcols = {'sub', 'trial', 'A', 'B', 'rt'};

%% Start simulation
s = 2; e = s + n - 1; % Counters
for repidx = 1:N
    tic
    %% Individual subject RTs are samples from the logNormal with parameter mu and sigma
    % Generate full dataset upfront
    simdata = []; % Table of stimulated data for each condition
    truevals = []; % Actual values for
    
    p = pset(mod(repidx-1, 5)+1);
    
    for sidx = 1:n
        %% Generate subject data
        if rand <= p; % Interaction occurs with probability p
            b = [beta, normrnd(interaction_mean, interaction_std)]; % Sample interaction from normal distribution
            int_present = 1;
        else
            b = [beta, normrnd(0, interaction_std)]; % Sample interaction from normal distribution
            int_present = 0;
        end
        
        % Arithmetic mean for RT (i.e., X * beta)
        m = simmeans(intercept, b);
        
        % Convert to the log normal parameters
        mu    = mufun(m, v);
        sigma = sigmafun(m, v);
        
        for cidx = 1:size(con, 1)
            f1 = con(cidx,1); % Factor 1 condition
            f2 = con(cidx,2); % Factor 2 condition
            rt = lognrnd(mu(cidx), sigma(cidx), M, 1); % Lognormal RTs
            simdata = [simdata; % Populate table of stimulated data for each condition
                sidx * ones(M, 1), (1:M)', f1 * ones(M, 1), f2 * ones(M, 1), rt];
        end
        truevals = [truevals, sidx, b, int_present];
    end
    groundtruth(repidx,:) = truevals; % factor1, factor2, and intercept values
    
    %% Now bootstrap data for individual and group analysis
    for bidx = 1:nBootstrap
        bootdata = [];
        for sidx = 1:n
            %% Estimate the beta weights from the lognormally distributed data for each subject using Weighted Least Squares
            subdata = [];
            for conIdx = 1:size(con, 1)
                sampdata = datasample(simdata(simdata(:,1) == sidx &...
                    simdata(:,3) == con(conIdx, 1) &...
                    simdata(:,4) == con(conIdx, 2), 5),...
                    nbootsamples, 'Replace', true); % Pull out subject data
                subdata = [subdata;...
                    sidx * ones(nbootsamples, 1),...
                    (1:nbootsamples)',...
                    con(conIdx, 1) * ones(nbootsamples, 1),...
                    con(conIdx, 2) * ones(nbootsamples, 1),...
                    sampdata];
            end
            bootdata = [bootdata; subdata];
        end
        
        %% Analyse individual bootstrap data
        tempout = [];
        for sidx = 1:n
            subdata = bootdata(bootdata(:,strcmp(subcols, 'sub')) == sidx, :);
            
            % Useful matrices
            % Intercept, factor 1, factor 2, interaction
            X = [ones(size(subdata, 1), 1), subdata(:,mstrfind(subcols, {'A', 'B'})), prod(subdata(:,mstrfind(subcols, {'A', 'B'})), 2)];
            y = subdata(:,strcmp(subcols, 'rt')); % y = X*Beta + error
            
            % Set up for using matlab fitglm
            dset = mat2dataset([y X(:,2:3)]); % Just include main effects
            dset.Properties.VarNames = {'y', 'A', 'B'};
            
            %% Weighted least squares
            modelspec = 'y ~ A*B'; % Full interaction model
            weights = 1./((X * (X\y)).^2); % Weights scale with 1/mean
            wls = fitglm(dset, modelspec, 'Distribution', 'normal', 'Link', 'identity', 'weights', weights);
            sigma_est_wls = std(wls.Residuals.Raw);
            
            %% Fit logNormal interaction model
            startparms = [wls.Coefficients.Estimate; sigma_est_wls]; % Start model at the weighted least squares estimate
            [fittedparms, fit, ~] = fminsearch(@(parms) lnL(parms, y, X), startparms, options); % Optimize using Nelder-Mead SIMPLEX algorithm
            beta_est_mle = fittedparms(1:4); % Beta weights
            s_est_mle    = fittedparms(5);   % Estimate of arithmetic standard deviation
            fittedMeans  = simmeans(fittedparms(1), fittedparms(2:4)); % Estimated of arithmetic mean
            
            %% Fit constrained model without interaction
            X_con = X(:,1:3); % Remove interaction
            startparms = [wls.Coefficients.Estimate(1:3); sigma_est_wls]; % Start model at the weighted least squares estimate
            [fittedparms_con, fit_con, ~] = fminsearch(@(parms) lnL_con(parms, y, X_con), startparms, options); % Optimize using Nelder-Mead SIMPLEX algorithm
            beta_est_mle_con = fittedparms(1:3); % Beta weights
            s_est_mle_con    = fittedparms(4);   % Estimate of arithmetic standard deviation
            fittedMeans_con  = simmeans_con(fittedparms_con(1), fittedparms_con(2:3)); % Estimated of arithmetic mean
            
            %% Compare full and constrained model
            g2   = -2 * (round(fit, 2)-round(fit_con, 2)); % G2
            pval = chi2pdf(g2, 1);                         % p-value
            
            % Store output
            tempout = [tempout, sidx, wls.Coefficients.Estimate(end), wls.Coefficients.pValue(end), beta_est_mle(end) pval];
        end
        bootout(bidx,:) = tempout;
        
        %% Group anlaysis on bootdata
        %% Now run a group-level ANOVA
        varnames = {'sub', 'A', 'B'};
        submeans = aggregate(bootdata, mstrfind(subcols, varnames), strcmp(subcols, 'rt'));  % Find average mean for each subject
        C = mat2cell(submeans(:,[1 2 3]), size(submeans, 1), ones(1, size(submeans(:,[1 2 3]), 2))); % Set up cell matrix for anovan
        [pv, table, stats, terms] = anovan(submeans(:,4), C, 'varnames', varnames, 'random', [1], 'model', 'full', 'display', 'off');
        
        anovaCols = table(1,:);
        interaction_line = table(strcmp(table(:,1), 'A*B'), :);
        groupout(bidx, :) = [interaction_line{strcmp(anovaCols, 'F')},...
            interaction_line{strcmp(anovaCols, 'd.f.')},...
            table{strcmp(table(:,1), 'sub*A*B'), strcmp(anovaCols, 'd.f.')},...
            interaction_line{strcmp(anovaCols, 'Prob>F')},...
            table{strcmp(table(:,1), 'sub*A*B'), strcmp(anovaCols, 'Mean Sq.')}];
        
        if mod(bidx, 100) == 0
            fprintf('Simulation %d\t', bidx);
            time = toc;
            displayTime(time);
            save(sprintf('largerNsimulation_%d', n))
        end
    end
    
    rboot = reshape(bootout', 5, n, size(bootout, 1));
    meanWlsEstimate = mean(rboot(2,:,:), 3);
    stdWlsEstimate  = std(rboot(2,:,:), [], 3);
    pCountWls       = mean(reshape(rboot(3,:,:), n, nBootstrap)' <= .05);
    meanFitEstimate = mean(rboot(4,:,:), 3);
    stdFitEstimate  = std(rboot(4,:,:), [], 3);
    pCountFit       = mean(reshape(rboot(5,:,:), n, nBootstrap)' <= .05);
    boSummary = [(1:n)', meanWlsEstimate', stdWlsEstimate', pCountWls', meanFitEstimate', stdFitEstimate', pCountFit'];
    
    bo(repidx,:) = reshape(boSummary', numel(boSummary), 1)';
    go(repidx,:) = mean(groupout(:,4) <= .05);
end
save(sprintf('largerNsimulation_%d', n))