%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Parameter Ranking via MC-Sobol %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% ***REQUIRES MATLAB STATISTICS TOOLBOX!!*** %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
tic;
format longG;
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 10000;

% (b) Initialize Parameter Values for Key Parameters with kp
%     The number of elements in kp should correspond to # of key parameters that are
%     considered in the CCU model. If the # of key params of the CCU model are 
%     changed, the length of kp should be adjusted accordingly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [1 1 1];                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              

% (c) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);

% (d) Define two identical parameter spaces for Sobol Analysis
%     If n_kp params are simulated n_sim times, a space of n_sim x n_kp is needed
ParamSpace = [];                     % Space for the Key Parameter of Interest
c_ParamSpace = [];                   % Space for Complementary Parameters to kp

% (e) Define Resolution of Kernel Distributions
np_1 = 10000;               % Resolution for Parameter KDs for MC sampling
np_2 = 50000;               % Resolution for Bayesian KDs for computing area difference


%% [2] Populate Parameter Space via MC Simulation
% Distributions can either be defined parametrically or non-parametrically. If 
% repeated experiment sampling under identical conditions is possible, then a  
% parametric sample distribution can be defined with mean and stdev. Otherwise, a PDF
% can be generated non-parametrically using sample data + Gaussian kernels. See below
% [EX] Parametric Distribution:
% ParamSpace(:,1) = 3.4.*randn(n_sim, 1) + 0.9                 >>> Mean=0.9, Stdev=3.4
% ParamSpace(:,3) = 10.*randn(n_sim, 1) + 24.1                 >>> Mean=24.1, Stdev=10
%
% [EX] Non-Parametric (Kernel) Distribution:
% data_var3 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2]..etc  >>> Sample data points
% [f_var3, xi_var3] = ksdensity(data_var3, 'npoints', np_1)    >>> Generate KDE PDF
% for i = 1:n_sim                                              >>> Sample from KDE PDF
%     ParamSpace(i,3) = randarb(xi_var3, f_var3);                  using "randarb"
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ParamSpace(:,1) = 25.*randn(n_sim, 1) + 95;                  % Parametric PDF for x1
c_ParamSpace(:,1) = 25.*randn(n_sim, 1) + 95;
data_x2 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2];
data_x3 = [82.4, 77.6, 83.3, 80.1];
[f_x2, xi_x2] = ksdensity(data_x2, 'npoints', np_1);
[f_x3, xi_x3] = ksdensity(data_x3, 'npoints', np_1);
for i = 1:n_sim                                               % Kernel PDFs for x2, x3
    % Need to sample twice to populate both ParamSpace and c_ParamSpace
    ParamSpace(i,2) = randarb(xi_x2, f_x2);
    c_ParamSpace(i,2) = randarb(xi_x2, f_x2);
    ParamSpace(i,3) = randarb(xi_x3, f_x3);
    c_ParamSpace(i,3) = randarb(xi_x3, f_x3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [3] Process Parameter Space for CCU Evaluation
% (a) Determine the # of CCU Model Outputs (# of Evaluation Metrics)
n_out = size(SampleModel(kp),2);

% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
outputs = [];               % Evaluated model outputs with MC sampled inputs
eval_P = [];                % kp of interest come from ParamSpace, 1st Sobol
eval_C = [];                % kp of interest come from ComplementarySpace, Total Sobol

% (c) Using Monte Carlo sampled parameter values in (2a), evaluate CCU Model
for i = 1:n_sim
    % Each Parameter Set is a Row in the ParamSpace matrix
    Parameter_Set = ParamSpace(i,:); 
    outputs(i,:) = SampleModel(Parameter_Set);
end

% (d) Compute Factors for 1st Order and Total Sobol
for i = 1:n_sim
    for j = 1:n_kp
        %%%% 1st Order Sobol Factors
        % Loop through all Parameters by indexing with j. For a certain parameter of
        % interest, (lets say 3rd parameter, j=3) kp becomes an n_sim x n_p matrix
        % where the column of values for the parameter of interest (3rd column) are
        % assigned from the Parameter Space. The columns for remaining parameters 
        % (j=/=3) are populated with values from the Complementary Space. 
        kp = [c_ParamSpace(i, 1:j-1), ParamSpace(i,j), c_ParamSpace(i, j+1:n_kp)];
        eval_P(i,:,j) = SampleModel(kp);
        %%%% Total Sobol Factors
        % Loop through all Parameters by indexing with j. For a certain parameter of
        % interest, (lets say 5th parameter, j=5) kp becomes an n_sim x n_p matrix
        % where the column of values for the parameter of interest (5th column) are
        % assigned from the Complementary Space. The columns for remaining parameters 
        % (j=/=5) are populated with values from the Parameter Space. 
        kp = [ParamSpace(i, 1:j-1), c_ParamSpace(i,j), ParamSpace(i, j+1:n_kp)];
        eval_C(i,:,j) = SampleModel(kp);
    end
end


%% [4] Computing Variances
% (a) Initialize Integrals and Total Variances
% Initialize Integral of Model Outputs via ANOVA decomposition (1 x n_out vector)
f0 = zeros(1, n_out);            
% Initialize Total Variance of Model Outputs initialized to zero (1 x n_out vector)
D = zeros(1, n_out);            

%(b) Estimating Total Output Variance, D, using Monte Carlo Simulated Model Outputs
for i = 1:n_out
    for j = 1:n_sim
    % The INTEGRAL OF THE MODEL OUTPUT is estimated by 1/N*SIGMA(f(x_k)) from k=1 to N
    % AKA, the AVERAGE of all Model Output values from the Monte Carlo simulations
    % EX: For j=1 (1st CCU model output), calculate the average output value
    f0(i) = f0(i) + outputs(j,i)/n_sim;
    % The TOTAL VARIANCE OF MODEL OUTPUT is estimated by SQUARING all outputs for a 
    % particular  
    % from Monte Carlo simulations then taking the average of SQUARED Model Outputs
    % Call this value D. The total variance is defined as D minus the square of the
    % average model outputs, aka, D - f0^2
    D(i) = D(i) + outputs(j,i).^2/n_sim;
    end
end

%(c) Estimating Total Output Variance, D, using Monte Carlo Simulated Model Outputs
for i = 1:n_out
    D(i) = D(i) - f0(i).^2;
end

%(d) Compute Partial Variances (1st Order Sobol)
% The Partial Variances for each Key Parameter j is defined as the Total Model Output 
% Variance D minus (1/2*simulation#) times the SQUARE of Monte Carlo Outputs with all 
% Parameters samples from Parameter Space minus the Monte Carlo Outputs with Parameter 
% j sampled from Parameter Space and non-j sampled from Complementary Space.
D_1st = [];
Factors_1st = 0;
for i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
        % D_1st is a n_kp x n_out vector. The values for Key Params are initialized 
        % to D(n_out) according to the CCU Model Outputs. For each MC simulation, the 
        % differences in output(k,j) and eval_P(k,j,i) are squared, then normalized by
        % dividing by [2*(n_sim)]. As we loop through n_sim, the Normalized Squared 
        % Difference for Key Param j is iteratively subtracted from the Total Model 
        % Output Variance D. When we finish subtracting the differences after looping 
        % through all n_sim, the Residual Variance is the Partial Variance
        Factors_1st = Factors_1st + ((outputs(k,j) - eval_P(k,j,i)).^2/(2*n_sim));
        end
        D_1st(i,j) = D(j) - Factors_1st;
        Factors_1st = 0;
    end
end

%(e) Compute Total Variances (Total Sobol)
D_Total = zeros(n_kp, n_out);
Factors_Total = 0;
for i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
        % D_Total is an n_kp x n_out vector that is initialized to zero. Unlike 1st
        % order Sobol, the parameter of interest is sampled from the Complementary
        % Space and not the Parameter Space. For each MC simulation, the differences
        % in output(k,j) and eval_C(k,j,i) are squared, then normalized by dividing by
        % [2*(n_sim)]. As we loop through n_sim, the Normalized Square Differences for
        % Key Param j is iteratively added. The total sum after iterative adding is
        % tha Total Variance contribution of Key Param j.
        Factors_Total = Factors_Total + ((outputs(k,j) - eval_C(k,j,i)).^2/(2*n_sim));
        end
        D_Total(i,j) = D_Total(i,j) + Factors_Total;
        Factors_Total = 0;
    end
end

%% [5] Determine Sobol Indices
% (a) Compute 1st Order and Total Sobol Indices
% For 1st Order Sobol, for each CCU model output, its D_1st is divided by D (total 
% output variance). The sum of 1st Order Sobol should equal 1 for the Conditional
% Independence Assumption to hold, however, should the sum (either 1st Order or Total)
% be >1, it means that parameter correlations are present. The greater the sum is to
% 1, the greater the parameter correlations and the less the conditional independence
% assumption holds
toc
Sobol_1st = [];
Sobol_Total = [];
for i = 1:n_out
    for j = 1:n_kp
    Sobol_1st(j,i) = D_1st(j,i)/D(i);
    Sobol_Total(j,i) = D_Total(j,i)/D(i);
    end
end

% (b) Rank the Parameters via Sobol Indices
for i = 1:n_out
    [score, rank] = sort(Sobol_1st(:,i), 'descend');
    fprintf('Order of Parameters by 1st Order Sobol Rank for CCU Eval Output %.4f \n', i)
    rank
    sprintf('The Rank Scores for the Above Order are: ')
    Sobol_1st(:,i)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [score2, rank2] = sort(Sobol_Total(:,i), 'descend');
    fprintf('Order of Parameters by Total Sobol Rank for CCU Eval Output %.4f \n', i)
    rank2
    sprintf('The Rank Scores for the Above Order are: ')
    Sobol_Total(:,i)
end

% (c) Export Results
xlswrite('MC-Sobol_ParameterSpace.xls', ParamSpace)
xlswrite('MC-Sobol_ComplementarySpace.xls', c_ParamSpace)
xlswrite('MC-Sobol_EvaluatedOutputs.xls', outputs)
xlswrite('MC-Sobol_TotVar+F0.xls', [D, f0])
xlswrite('MC-Sobol_Variance_1st.xls', D_1st)
xlswrite('MC-Sobol_Variance_Total.xls', D_Total)
xlswrite('MC-Sobol_1stOrderSobol.xls', [score, rank])
xlswrite('MC-Sobol_TotalSobol.xls', [score2, rank2])
