%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PARAMETER RANKING via MC-SOBOL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
tic;
format longG;
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 1000;

% (b) Initialize Parameter Values for Key Parameters with kp
%     The number of elements in kp should correspond to # of key parameters 
%     that are considered in the model. If the # of key params of the model 
%     are changed, the length of kp should be adjusted accordingly.
%     GUI Support will be added in the future
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kp = [1 1 1 1 1 1 1 1 1 1];                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

% (c) Scalar # of Key Parameters. Count columns of kp array.
n_kp = size(kp, 2);

% (d) Define two identical parameter spaces for Sobol Analysis
%     If n_kp params are simulated n_sim times, space is n_sim x n_kp 
ParamSpace = [];                 % Space for the Key Parameter of Interest
c_ParamSpace = [];               % Space for Complementary Parameters to kp



%% [2] Populate Parameter Space via MC Simulation
% Distributions can either be defined parametrically or non-parametrically. 
% If  repeated experiment sampling under identical conditions is possible, 
% then a parametric sample distribution can be defined with mean and stdev.
% Otherwise, a PDF can be generated non-parametrically via kernel density
% [EX] Parametric Distribution:
% ParamSpace(:,1) = 3.4.*randn(n_sim, 1) + 0.9                 
% ParamSpace(:,3) = 10.*randn(n_sim, 1) + 24.1                
%
% [EX] Kernel Distribution using fitdist
% data_x1 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2]  % Sample data
%

% [EX] Kernel Distribution using ksdensity (slower than fitdist!)
% ParamSpace(:,1) = 35.*randn(n_sim, 1) + 325;      % Parametric PDF for x1
% data_x2 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2];
% data_x3 = [82.4, 77.6, 83.3, 80.1];
% [f_x2, xi_x2] = ksdensity(data_x2, 'npoints', np_1);
% [f_x3, xi_x3] = ksdensity(data_x3, 'npoints', np_1);
% for i = 1:n_sim                                               % Kernel PDFs for x2, x3
%     ParamSpace(i,2) = randarb(xi_x2, f_x2);
%     ParamSpace(i,3) = randarb(xi_x3, f_x3);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
ParamSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;                        % Parametric PDF for avg_GPC
ParamSpace(:,8) = 13030.6.*randn(n_sim, 1) + 65153;                        % Parametric PDF for C_PBR
c_ParamSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;                      % Parametric PDF for avg_GPC
c_ParamSpace(:,8) = 13030.6.*randn(n_sim, 1) + 65153;                      % Parametric PDF for C_PBR

% [List Kernel Data Points]
data_x2 = [0.89, 1.07, 1.24, 1.00, 0.4, 0.8, 0.84, 0.81, 1.11, 0.82, 1];   % Kernel Data for Mu_max
data_x3 = [3.6, 2.042, 1.4];                                               % Kernel Data for kLa
data_x4 = [0.3, 0.14, 0.18, 0.22, 0.1218, 0.19, 0.201, 0.225];      % Kernel Data for y_Lipid
data_x5 = [0.1358, 0.1184, 0.1283, 0.116];                                 % Kernel Data for P_CO2
data_x6 = [0.402, 0.424, 0.465, 0.515, 0.562, 0.592, 0.599, 0.581, 0.542, 0.493, 0.445, 0.411]; % Kernel Data for Daylight %
data_x7 = [0.95, 0.93, 0.99, 0.97];                                        % Kernel Data for Har+Dew Recovery
data_x9 = [187.5, 247.2, 278.06, 313.88, 469.17, 313.06, 288.88, 250, 266.11, 241.66]; % Kernel Data for Coal Grid Elec. GWI
data_x10 = [3.234, 4.62, 6.006, 3.469, 3.528, 3.527, 6.2, 2.74, 1.59, 1.13];           % Kernel Data for N-Fertilizer GWI

% [Define Kernel Distributions]
f_x2 = fitdist(data_x2', 'Kernel', 'Support', [0.01, 10]);
f_x3 = fitdist(data_x3', 'Kernel', 'Support', [0.01, 100]);
f_x4 = fitdist(data_x4', 'Kernel', 'Support', [0.01, 0.9]);
f_x5 = fitdist(data_x5', 'Kernel', 'Support', [0.01, 0.9]);
f_x6 = fitdist(data_x6', 'Kernel', 'Support', [0.25, 0.75]);
f_x7 = fitdist(data_x7', 'Kernel', 'Support', [0.7, 1]);
f_x9 = fitdist(data_x9', 'Kernel', 'Support', [1, 999]);
f_x10 = fitdist(data_x10', 'Kernel', 'Support', [0.01, 100]);

% [Sample from Kernel Distributions]
% Random Sample to Populate Parameter Space
ParamSpace(:,2) = random(f_x2, n_sim, 1);
ParamSpace(:,3) = random(f_x3, n_sim, 1);
ParamSpace(:,4) = random(f_x4, n_sim, 1);
ParamSpace(:,5) = random(f_x5, n_sim, 1);
ParamSpace(:,6) = random(f_x6, n_sim, 1);
ParamSpace(:,7) = random(f_x7, n_sim, 1);
ParamSpace(:,9) = random(f_x9, n_sim, 1);
ParamSpace(:,10) = random(f_x10, n_sim, 1);
% Random Sample to Populate Complementary Parameter Space
c_ParamSpace(:,2) = random(f_x2, n_sim, 1);
c_ParamSpace(:,3) = random(f_x3, n_sim, 1);
c_ParamSpace(:,4) = random(f_x4, n_sim, 1);
c_ParamSpace(:,5) = random(f_x5, n_sim, 1);
c_ParamSpace(:,6) = random(f_x6, n_sim, 1);
c_ParamSpace(:,7) = random(f_x7, n_sim, 1);
c_ParamSpace(:,9) = random(f_x9, n_sim, 1);
c_ParamSpace(:,10) = random(f_x10, n_sim, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% [3] Process Parameter Space for CCU Evaluation
% (a) Determine the # of CCU Model Outputs (# of Evaluation Metrics)
n_out = size(CCUS_Biocrude(kp),2);

% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
outputs = zeros(n_sim, n_out);   % Evaluated outputs with MC sampled inputs
eval_P = zeros(n_sim, n_out);    % kp of interest come from ParamSpace, 1st Sobol
eval_C = zeros(n_sim, n_out);    % kp of interest come from ComplementarySpace, Total Sobol

% (c) Using Monte Carlo sampled parameter values in (2a), evaluate CCU Model
for i = 1:n_sim
    % Each Parameter Set is a Row in the ParamSpace matrix
    Parameter_Set = ParamSpace(i,:); 
    outputs(i,:) = CCUS_Biocrude(Parameter_Set);
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
        eval_P(i,:,j) = CCUS_Biocrude(kp);
        %%%% Total Sobol Factors
        % Loop through all Parameters by indexing with j. For a certain parameter of
        % interest, (lets say 5th parameter, j=5) kp becomes an n_sim x n_p matrix
        % where the column of values for the parameter of interest (5th column) are
        % assigned from the Complementary Space. The columns for remaining parameters 
        % (j=/=5) are populated with values from the Parameter Space. 
        kp = [ParamSpace(i, 1:j-1), c_ParamSpace(i,j), ParamSpace(i, j+1:n_kp)];
        eval_C(i,:,j) = CCUS_Biocrude(kp);
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
D_1st = zeros(n_kp, n_out);
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
Sobol_1st = zeros(n_kp, n_out);
Sobol_Total = zeros(n_kp, n_out);
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
toc
