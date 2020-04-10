%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% PARAMETER RANKING via MC-SOBOL %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Jeehwan Lee, KAIST, stevelee@kaist.ac.kr %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clc;
tic;
format longG;
progress = waitbar(0, '1', 'Name', 'Running PR-MC-Sobol...');
total_steps = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 10000;

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
ParSpace = [];                 % Space for the Key Parameter of Interest
c_ParSpace = [];               % Space for Complementary Parameters to kp

%%Progress Bar%%
waitbar(20/total_steps, progress, 'Generating Parameter Spaces...');


%% [2] Populate Parameter Space via MC Simulation
% Distributions can either be defined parametrically or non-parametrically. 
% If  repeated experiment sampling under identical conditions is possible, 
% then a parametric sample distribution can be defined with mean and stdev.
% Otherwise, a PDF can be generated non-parametrically via kernel density
% [EX] Parametric Distribution:
% ParSpace(:,1) = 3.4.*randn(n_sim, 1) + 0.9                 
% ParSpace(:,3) = 10.*randn(n_sim, 1) + 24.1                
%
% [EX] Kernel Distribution using fitdist
% data_x1 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2]  % Sample data
%
% [EX] Kernel Distribution using ksdensity (slower than fitdist!)
% ParSpace(:,1) = 35.*randn(n_sim, 1) + 325;     % Parametric PDF for x1
% data_x2 = [13.5, 16.4, 14.4, 19.6, 15.0, 15.9, 16.2];
% data_x3 = [82.4, 77.6, 83.3, 80.1];
% [f_x2, xi_x2] = ksdensity(data_x2, 'npoints', np_1);
% [f_x3, xi_x3] = ksdensity(data_x3, 'npoints', np_1);
% for i = 1:n_sim                                  % Kernel PDFs for x2,x3
%     ParSpace(i,2) = randarb(xi_x2, f_x2);
%     ParSpace(i,3) = randarb(xi_x3, f_x3);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Define Parametric Distributions for Key Params]
ParSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;               % avg_GPC
ParSpace(:,8) = 13030.6.*randn(n_sim, 1) + 65153;               % C_PBR
c_ParSpace(:,1) = 0.0016.*randn(n_sim, 1) + 0.0224;             % avg_GPC
c_ParSpace(:,8) = 13030.6.*randn(n_sim, 1) + 65153;             % C_PBR

% [List Kernel Data Points]
data_x2 = [0.89, 1.07, 1.24, 1.00, 0.4, 0.8, 0.84,...
           0.81, 1.11, 0.82, 1];           % Kernel Data for Mu_max
data_x3 = [3.6, 2.042, 1.4];               % Kernel Data for kLa
data_x4 = [0.3, 0.14, 0.18, 0.22, 0.1218,...
           0.19, 0.201, 0.225];            % Kernel Data for y_Lipid
data_x5 = [0.1358, 0.1184, 0.1283, 0.116]; % Kernel Data for P_CO2
data_x6 = [0.402, 0.424, 0.465, 0.515, 0.562, 0.592, 0.599, 0.581,...
           0.542, 0.493, 0.445, 0.411];    % Kernel Data for Daylight %
data_x7 = [0.95, 0.93, 0.99, 0.97];        % Kernel Data for Recovery %
data_x9 = [187.5, 247.2, 278.06, 313.88, 469.17, 313.06, 288.88, 250,...
           266.11, 241.66];                % Kernel Data for Elec. GWI
data_x10 = [3.234, 4.62, 6.006, 3.469, 3.528, 3.527, 6.2, 2.74, 1.59,...
            1.13];                         % Kernel Data for N-Fert. GWI

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
ParSpace(:,2) = random(f_x2, n_sim, 1);
ParSpace(:,3) = random(f_x3, n_sim, 1);
ParSpace(:,4) = random(f_x4, n_sim, 1);
ParSpace(:,5) = random(f_x5, n_sim, 1);
ParSpace(:,6) = random(f_x6, n_sim, 1);
ParSpace(:,7) = random(f_x7, n_sim, 1);
ParSpace(:,9) = random(f_x9, n_sim, 1);
ParSpace(:,10) = random(f_x10, n_sim, 1);
% Random Sample to Populate Complementary Parameter Space
c_ParSpace(:,2) = random(f_x2, n_sim, 1);
c_ParSpace(:,3) = random(f_x3, n_sim, 1);
c_ParSpace(:,4) = random(f_x4, n_sim, 1);
c_ParSpace(:,5) = random(f_x5, n_sim, 1);
c_ParSpace(:,6) = random(f_x6, n_sim, 1);
c_ParSpace(:,7) = random(f_x7, n_sim, 1);
c_ParSpace(:,9) = random(f_x9, n_sim, 1);
c_ParSpace(:,10) = random(f_x10, n_sim, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Progress Bar%%
waitbar(100/total_steps, progress, 'Initializing Monte Carlo Simulations');


%% [3] Process Parameter Space for CCU Evaluation
% (a) Determine the # of CCU Model Outputs (# of Evaluation Metrics)
n_out = size(CCUS_Biocrude(kp),2);

% (b) Create containers for Evaluated Model Outputs
% NOTE: eval_P and eval_C are n_sim x n_out x n_kp matrix. 
outputs = zeros(n_sim, n_out);   % Evaluated outputs with MC sampled inputs
eval_P = zeros(n_sim, n_out);    % kp from ParSpace, 1st Sobol
eval_C = zeros(n_sim, n_out);    % kp from ComplementarySpace, Total Sobol

% (c) Using Monte Carlo sampled parameter values in (2a), evaluate model
for i = 1:n_sim
    % Each Parameter Set is a Row in the ParSpace matrix
    Parameter_Set = ParSpace(i,:); 
    outputs(i,:) = CCUS_Biocrude(Parameter_Set);
end

% (d) Compute Factors for 1st Order and Total Sobol
for i = 1:n_sim
    for j = 1:n_kp
        %%%% 1st Order Sobol Factors
        % Loop through all Parameters by indexing with j. For a certain 
        % parameter of interest, (lets say 3rd parameter, j=3) kp becomes 
        % a n_sim x n_p matrix where the column values for the parameter of
        % interest (3rd column) are assigned from the Parameter Space. The 
        % columns for remaining parameters (j=/=3) are populated with 
        % values from the Complementary Space. 
        kp = [c_ParSpace(i,1:j-1), ParSpace(i,j), c_ParSpace(i,j+1:n_kp)];
        eval_P(i,:,j) = CCUS_Biocrude(kp);
        %%%% Total Sobol Factors
        % Loop through all Parameters by indexing with j. For a certain 
        % parameter of interest, (lets say 5th parameter, j=5) kp becomes a
        % n_sim x n_p matrix where the column values for the parameter of 
        % interest (5th column) are assigned from the Complementary Space. 
        % The columns for remaining parameters (j=/=5) are populated with 
        % values from the Parameter Space. 
        kp = [ParSpace(i,1:j-1), c_ParSpace(i,j), ParSpace(i,j+1:n_kp)];
        eval_C(i,:,j) = CCUS_Biocrude(kp);
    end
    %%Progress Bar%%
    waitbar((100 + ((i/n_sim)*820))/total_steps, progress,...
             sprintf('Processing %d of %d Simulations', i, n_sim));
end
%%Progress Bar%%
waitbar(920/total_steps, progress, 'Computing Sobol Variances');


%% [4] Computing Variances
% (a) Initialize Integrals and Total Variances
% Initialize Integral of Model Outputs via ANOVA decomposition 
f0 = zeros(1, n_out);            
% Initialize Total Variance of Model Outputs initialized to zero 
D = zeros(1, n_out);            

% (b) Estimating Total Output Variance, D, using MC Model Outputs
for i = 1:n_out
    for j = 1:n_sim
        % The Integral of Model Outputs is estimated by 1/N*SIGMA(f(x_k)) 
        % from k=1 to N, AKA, the average of all Model Output values from 
        % the Monte Carlo simulations
        % EX: For j=1 (1st CCU model output), calculate the average output 
        % value f0(i) = f0(i) + outputs(j,i)/n_sim;
        % The TOTAL VARIANCE of Model Outputs is estimated by squaring all 
        % outputs for an MC simulation then taking the average of squared
        % Model Outputs. Call this value D. The total variance is defined 
        % as D minus the square of the average model outputs, D - f0^2
        D(i) = D(i) + outputs(j,i).^2/n_sim;
    end
end

% (c) Estimating Total Output Variance, D, using MC Model Outputs
for i = 1:n_out
    D(i) = D(i) - f0(i).^2;
end

% (d) Compute Partial Variances (1st Order Sobol)
% The Partial Variances for each Uncertain Parameter j is defined as the 
% Total Model Output Variance D minus (1/2*simulation#) times the square 
% of Monte Carlo Outputs with all parameters samples from Parameter Space 
% minus the Monte Carlo Outputs with Parameter j sampled from Parameter 
% Space and non-j sampled from Complementary Space.
D_1st = zeros(n_kp, n_out);
F_1st = 0;                                 % Initialize 1st Order Factors
for i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
        % D_1st is a n_kp x n_out vector. The values for Key Params are 
        % initialized to D(n_out) according to the CCU Model Outputs. For 
        % each MC simulation, output(k,j)-eval_P(k,j,i) is squared, then 
        % normalized by dividing by [2*(n_sim)]. As we loop through n_sim, 
        % the Normalized Squared Difference for Key Param j is iteratively 
        % subtracted from the Total Model Output Variance D. When we finish 
        % subtracting the differences after looping through all n_sim, the 
        % Residual Variance is the Partial Variance for that parameter
        F_1st = F_1st + ((outputs(k,j) - eval_P(k,j,i)).^2/(2*n_sim));
        end
        D_1st(i,j) = D(j) - F_1st;
        F_1st = 0;
    end
end

%(e) Compute Total Variances (Total Sobol)
D_Total = zeros(n_kp, n_out);
F_Total = 0;                               % Initialize Total Sobol Factors
for i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
        % D_Total is an n_kp x n_out vector that is initialized to zero. 
        % Unlike 1st order Sobol, the parameter of interest is sampled from 
        % the Complementary Space. For each MC simulation, the differences
        % in output(k,j) and eval_C(k,j,i) are squared, then normalized by 
        % dividing by [2*(n_sim)]. As we loop through n_sim, the Normalized 
        % Square Differences for KeyParam j is iteratively added. The total 
        % sum after iterative adding is the Total Variance contribution of 
        % Key Param j.
        F_Total = F_Total + ((outputs(k,j) - eval_C(k,j,i)).^2/(2*n_sim));
        end
        D_Total(i,j) = D_Total(i,j) + F_Total;
        F_Total = 0;
    end
end
%%Progress Bar%%
waitbar(950/total_steps, progress, 'Exporting Sobol Indices');


%% [5] Determine Sobol Indices
% (a) Compute 1st Order and Total Sobol Indices
% For 1st Order Sobol, for each CCU model output, its D_1st is divided by D
% (total output variance). The sum of 1st Order Sobol should equal 1 for 
% the Conditional Independence Assumption to hold, however, should the sum 
% (either 1st Order or Total) be >1, it means that parameter correlations 
% are present. The greater the sum is to 1, the greater the parameter 
% correlations and the less the conditional independence assumption holds
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
    % 1st Order Indices
    [score, rank] = sort(Sobol_1st(:,i), 'descend');
    fprintf('1st Order Sobol Ranks for CCU Eval Output %.4f \n', i)
    rank
    sprintf('The Rank Scores for the Above Order are: ')
    Sobol_1st(:,i)
    fprintf('Sum of 1st Order Sobol Indices for Output %.4f \n', i)
    sum(Sobol_1st(:,i))
    % Total Sobol Indices
    [score2, rank2] = sort(Sobol_Total(:,i), 'descend');
    fprintf('Total Sobol Ranks for CCU Eval Output %.4f \n', i)
    rank2
    sprintf('The Rank Scores for the Above Order are: ')
    Sobol_Total(:,i)
    fprintf('Sum of Total Sobol Indices for CCU Eval Output %.4f \n', i)
    sum(Sobol_Total(:,i))
end


%% [6] Plot and Export Results
% (a) Plot Distributions of Model Outputs
nbins = 30;             % Define resolution of output metric's historgram

% (b) Generate plots (Example for via For loop below)
%for i = 1:n_out
    %figure(i)
    %histogram(outputs(:,i), nbins, 'facecolor', [1/i, 1/(4*i), 1/(16*i)])
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Output #1
    figure('Name', 'Variance of Unit Prod Cost')
    histogram(outputs(:,1), nbins, 'facecolor', [0, 0, 0])
    xlabel('Biocrude Unit Prod. Cost, $/kg')
    ylabel('Frequency')
    % Output #2
    figure('Name', 'Variance of Specific GWI')
    histogram(outputs(:,2), nbins, 'facecolor', [0, 0.5, 0])
    xlabel('Specific GWI, kg-CO2-eq/kg')
    ylabel('Frequency')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (c) Export Results
xlswrite('MC-Sobol_ParameterSpace.xls', ParSpace)
xlswrite('MC-Sobol_ComplementarySpace.xls', c_ParSpace)
xlswrite('MC-Sobol_EvaluatedOutputs.xls', outputs)
xlswrite('MC-Sobol_TotVar+F0.xls', [D, f0])
xlswrite('MC-Sobol_Variance_1st.xls', D_1st)
xlswrite('MC-Sobol_Variance_Total.xls', D_Total)
xlswrite('MC-Sobol_1stOrderSobol.xls', [score, rank])
xlswrite('MC-Sobol_TotalSobol.xls', [score2, rank2])
%%Progress Bar%%
waitbar(1000/total_steps, progress, 'Complete');
delete(progress)
toc