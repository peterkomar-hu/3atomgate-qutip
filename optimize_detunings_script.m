% This script finds the optimal values of Delta1 and Delta2, which minimize
% the DEVIATION function (from "deviation.m" in the same directory)
% It does so for a list of C values, and given kappa, gamma, and beta
% parameters
%
% C = cooperativity of the atom-cavity system
% kappa = cavity decay rate
% gamma = atomic excited state decay rate
% beta = Omega / (gamma * sqrt(C)), where Omega is the driving laser's strength
% Delta1 = Delta_E = detuning of the control atom's excited state
% Delta2 = Delta_e = detuning of the qubit atoms' excited state

clear all;

% parameters
kappa = 1;  % cavity dacay rate
gamma = 0.01;   % atomic excited state decay rate
C_list = [5, 50, 500, 5000]; % list of cooperativities to run the optimization
beta_list = [0.25, 0.125];   % beta: ratio between the driving strength and sqrt(C) * gamma

% output containers
Delta_E_opt = zeros(length(C_list), length(beta_list) );    % optimal values for Delta1
Delta_e_opt = zeros(length(C_list), length(beta_list) );    % optimal values for Delta2
K_opt = zeros(length(C_list), length(beta_list) );          % minimal value of deviatio

% open output files
file1 = fopen('Delta1.txt', 'w');
file2 = fopen('Delta2.txt', 'w');
fileK = fopen('K_min.txt', 'w');

% write headers to output files: the list of beta values
for j = 1:length(beta_list)
    fprintf( file1, '\t%f', beta_list(j) );
    fprintf( file2, '\t%f', beta_list(j) );
    fprintf( fileK, '\t%f', beta_list(j) );
end
fprintf( file1, '\n');
fprintf( file2, '\n');
fprintf( fileK, '\n');

% main loop, running through the list of C values
for i =  1:length(C_list)
    
    % print the value of C as the first column of all three files
    fprintf( file1, '%f', C_list(i) );
    fprintf( file2, '%f', C_list(i) );
    fprintf( fileK, '%f', C_list(i) );
    
    % for each beta value
    for j = 1:length(beta_list)
        
        C = C_list(i);
        
        % strength of the drive
        Omega = gamma * sqrt( C ) * beta_list(j);  
        
        % theoretical optima for the detunings (Delta1 and Delta2)
        Delta_E_th = gamma/2 * sqrt(4*C + 1);
        Delta_e_th = C * gamma / sqrt(4 * C + 1);

        % find minimum of DEVIATION, by finding the optimal values of Delta1 and
        % Delta2, using FMINSEARCH
        K = @(d) deviation(C, Omega, kappa, gamma, d(1), d(2));
        [Delta_opt, K_min, exitflag] = ...
            fminsearch(K, [Delta_E_th, Delta_e_th], optimset('MaxIter',1e7)); 
        
        % checking for convergence
        if exitflag ~= 1
            error_msg = ...
            sprintf('fminsearch did not converge for C = %f, Omega = %f',...
                C, Omega)
        end
        
        % save the outputs
        Delta_E_opt(i,j) = Delta_opt(1);
        Delta_e_opt(i,j) = Delta_opt(2);
        K_opt(i,j) = K_min;
        
        % print the outputs to file
        fprintf( file1, '\t%.16e', Delta_opt(1) );
        fprintf( file2, '\t%.16e', Delta_opt(2) );
        fprintf( fileK, '\t%.16e', K_min );
        
    end
    
    % final end-of-line
    fprintf( file1, '\n' );
    fprintf( file2, '\n' );
    fprintf( fileK, '\n' );
    
end

% clos files
fclose(file1);
fclose(file2);
fclose(fileK);












