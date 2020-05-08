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
progress = waitbar(0, 'Running...', 'Name', 'Running PR-MC-Sobol...');
set(findall(progress), 'Units', 'Normalized')
set(progress, 'Position', [0.25, 0.4, 0.18, 0.12])
total_steps = 1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% [1] Initialize Parameters and Spaces
% (a) Define # of Model Evaluation (# of MC sampling of Parameter Sets)
n_sim = 5000;
export = 1;
if n_sim > 50000
    simwarn = questdlg('WARNING: For n_sim >50,000 there could be issues with exporting results via xlswrite. Disable xlswrite?',...
        'WARNING',...
        'No','Yes','Yes');
    switch simwarn
        case 'Yes'
            export = 0;
        case 'No'
            export = 1;
    end
end
            
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
fx = zeros(n_sim, n_out);  % Evaluated outputs with all Params from ParSpace
fx_P = zeros(n_sim, n_out);% Evaluated outputs with i from ParSpace, ~i from c_ParSpace
fx_C = zeros(n_sim, n_out);% Evaluated outputs with i from c_ParSpacek, ~i from ParSpace

% (c) Evaluate Model from Monte Carlo Sampled Inputs (ParSpace)
for i = 1:n_sim
    % Each Parameter Set is a Row in the ParSpace matrix
    Parameter_Set = ParSpace(i,:); 
    fx(i,:) = CCUS_Biocrude(Parameter_Set);
end

% (d) Generate Function Output Space based on i and ~i
for i = 1:n_sim
    for j = 1:n_kp
        %%%% fx_P = f(x_ik, x'_~ik)
        % Loop through all Parameters by indexing with j. For a certain 
        % parameter of interest, (lets say 3rd parameter, j=3) kp becomes 
        % a n_sim x n_p matrix where the column values for the parameter of
        % interest (3rd column) are assigned from the Parameter Space. The 
        % columns for remaining parameters (j=/=3) are populated with 
        % values from the Complementary Space. This input parameter set
        % (x_ik, x'_~ik) is evaluated to compute D_i (partial variance for 
        % parameter i)
        kp = [c_ParSpace(i,1:j-1), ParSpace(i,j), c_ParSpace(i,j+1:n_kp)];
        fx_P(i,:,j) = CCUS_Biocrude(kp);
        %%%% fx_C = f(x'_ik, x_~ik)
        % Loop through all Parameters by indexing with j. For a certain 
        % parameter of interest, (lets say 5th parameter, j=5) kp becomes a
        % n_sim x n_p matrix where the column values for the parameter of 
        % interest (5th column) are assigned from the Complementary Space. 
        % The columns for remaining parameters (j=/=5) are populated with 
        % values from the Parameter Space. This input parameter set
        % (x'_ik, x_~ik) is evaluated to compute D_Tot_i (total variance 
        % involving parameter i)
        kp = [ParSpace(i,1:j-1), c_ParSpace(i,j), ParSpace(i,j+1:n_kp)];
        fx_C(i,:,j) = CCUS_Biocrude(kp);
    end
    %%Progress Bar%%
    waitbar((100 + ((i/n_sim)*820))/total_steps, progress,...
    sprintf('Processing %d of %d Simulations', i, n_sim));
end
%%Progress Bar%%
waitbar(920/total_steps, progress, 'Computing Sobol Variances');



%% [4] Computing Variances
% (a) Initialize Integrals and Total Variances
% Initialize the Integral of Model Outputs (f0^2) 
f0 = zeros(1, n_out);            
% Initialize Total Variance of Model Outputs 
D = zeros(1, n_out);            

% (b) Compute the Average of Model Outputs, f0
for i = 1:n_out 
    for j = 1:n_sim
        % The average of the Model Outputs is defined as the sum of all
        % model outputs for each output in fx, divided by the number of
        % Monte Carlo Simulations
        f0(i) = f0(i) + fx(j,i);
    end
    f0(i) = f0(i)/n_sim;
end

% (c) Estimating Total Output Variance, D, using MC Model Outputs
for i = 1:n_out
    for j = 1:n_sim
        % The TOTAL variance of Model Outputs is defined using ANOVA, and
        % is the integral of the square of function outputs minus the
        % square of the average of the function outputs. This integral is
        % estimated numerically via Monte Carlo as follows:
        D(i) = D(i) + ((fx(j,i)^2)-(f0(i)^2));
    end
    D(i) = D(i)/n_sim;
end

% (d) Compute Partial Variances (1st Order Sobol)
% The Partial Variances for each Uncertain Parameter j is defined as the 
% total variance D minus (1/2N) times the square of (fx - fx_P)
D_1st = zeros(n_kp, n_out);
F_Par = zeros(n_kp, n_out);           % Initialize Partial Variance Factors
for i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
            F_Par(i,j) = F_Par(i,j) + (fx(k,j) - fx_P(k,j,i))^2;
        end
        D_1st(i,j) = D(j) - (F_Par(i,j)/(n_sim*2));
    end
end

%(e) Compute Total Variances (Total Sobol)
D_Tot = zeros(n_kp, n_out);
F_cPar = zeros(n_kp, n_out);          % Initialize Total Sobol Factors
for i = 1:n_kp
    for j = 1:n_out
        for k = 1:n_sim
            % Unlike 1st order Sobol, the parameter of interest is sampled
            % from the Complementary Space. For each MC simulation, fx(k,j)
            % - fx_C(k,j,i) is squared then normalized by dividing by 
            % 2(n_sim). As we loop through n_sim, the Normalized Square 
            % Differences for KeyParam jis iteratively added
            F_cPar(i,j) = F_cPar(i,j) + (fx(k,j) - fx_C(k,j,i))^2;
        end
        D_Tot(i,j) = F_cPar(i,j)/(n_sim*2);
    end
end
%%Progress Bar%%
waitbar(950/total_steps, progress, 'Exporting Sobol Indices');



%% [5] Determine Sobol Indices
% (a) Compute Sum of Partial and Total Variances
% The Sum of Partial Variances is computed by adding D_1st/D for each
% particular output j. The Sum of Total Variances is computed by adding
% D_Tot/D for each output j. Note, if the Naive assumption for Sobol Method
% truly holds, the sum of all D_1st/D should not exceed 1. However, if
% the model parameters are positively/negatively correlated, the sum of 
% partial and total indices can be >1 and <1, respectively.
Sum_Partial = zeros(n_out);
Sum_Total = zeros(n_out);
for i = 1:n_out
    Sum_Partial(i) = sum(D_1st(:,i));
    Sum_Total(i) = sum(D_Tot(:,i));
end

% (b) Compute Rank Scores based on Partial/Total Variances
% The closer the sum of partial variances is to the total model output
% variances, the more additive the model is. 
Sobol_1st = zeros(n_kp, n_out);
Sobol_Total = zeros(n_kp, n_out);
for i = 1:n_out
    for j = 1:n_kp
        Sobol_1st(j,i) = D_1st(j,i)/D(i);
        Sobol_Total(j,i) = D_Tot(j,i)/D(i);
    end
end

% (b) Rank the Parameters via Sobol Indices
for i = 1:n_out
    % 1st Order Sobol Indices
    [score, rank] = sort(Sobol_1st(:,i), 'descend');
    fprintf('1st Order Sobol Ranks for CCU Eval Output %.0f \n', i)
    rank
    sprintf('The Partial Variance Indices (1st Order Sobol) for the Above Order are: ')
    Sobol_1st(:,i)
    fprintf('Sum of Partial Variance Indices for Output %.0f \n', i)
    sum(Sobol_1st(:,i))
    fprintf('Relative 1st Order Sobol Indices for Output %.0f \n', i)
    Sobol_1st(:,i)*D(i)/Sum_Partial(i)
    % Total Sobol Indices
    [score2, rank2] = sort(Sobol_Total(:,i), 'descend');
    fprintf('Total Sobol Ranks for CCU Eval Output %.0f \n', i)
    rank2
    sprintf('The Total Variance Indices (Total Sobol) for the Above Order are: ')
    Sobol_Total(:,i)
    fprintf('Sum of Total Sobol Indices for CCU Eval Output %.0f \n', i)
    sum(Sobol_Total(:,i))
    fprintf('Relative Total Sobol Indices for Output %.0f \n', i)
    Sobol_Total(:,i)*D(i)/Sum_Total(i)
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
    histogram(fx(:,1), nbins, 'facecolor', [0, 0, 0])
    xlabel('Biocrude Unit Prod. Cost, $/kg')
    ylabel('Frequency')
    % Output #2
    figure('Name', 'Variance of Specific GWI')
    histogram(fx(:,2), nbins, 'facecolor', [0, 0.5, 0])
    xlabel('Specific GWI, kg-CO2-eq/kg')
    ylabel('Frequency')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (c) Export Results
%%Progress Bar%%
waitbar(1000/total_steps, progress, 'Complete');
delete(progress)
toc

if export == 1
    xlswrite('MC-Sobol_EvaluatedOutputs.xls', fx)
    xlswrite('MC-Sobol_TotVar+F0.xls', [D, f0])
    xlswrite('MC-Sobol_Variance_1st.xls', D_1st)
    xlswrite('MC-Sobol_Variance_Total.xls', D_Tot)
    xlswrite('MC-Sobol_1stOrderSobol.xls', [score, rank])
    xlswrite('MC-Sobol_TotalSobol.xls', [score2, rank2])
end