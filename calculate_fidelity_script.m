% This script calculates the the fidelities for each density
% matrix rho in time series, read from tab-delimited files, named
% "data_rho_####.txt", and writes the fidelity time series to "fidelity_####.txt".
% 
% Fideltiy is defeined as the abs square of the maximal overlap between the actual
% state and 0.5 * ( |00> + e^(i*phi2)*|01> + e^(i*phi1)*|10> - e^(i*(phi1 + phi2))|11> ).
% as phi1 and phi2 are optimized to give the maximal overlap


clear all;

% the list of registry numbers to specify the input files
reglist = linspace(3221,3224,4);

% for all registry numbers
for m = 1:length(reglist)
    reg = reglist(m)
    
    % load the entire file as a matrix
    regstr = num2str(reg,'%04.0f');
    rhodata = importdata(['./results/data_rho_',regstr,'.txt']);
    
    % read timepoint from first column
    time = rhodata(:,1);
    
    % read entries of the density matrix:
    %   (diagonal elements are read from their respective position, as from a
    %   flattened matrix,
    %   real (imag) parts of the off-digonal elements are read from the upper
    %   (lower) triangle of the (flattened) matrix)
    rho00 = rhodata(:,2);
    rho01 = rhodata(:,3) + 1j* rhodata(:,6);
    rho02 = rhodata(:,4) + 1j* rhodata(:,10);
    rho03 = rhodata(:,5) + 1j* rhodata(:,14);

    rho10 = rhodata(:,3) - 1j* rhodata(:,6);
    rho11 = rhodata(:,7);
    rho12 = rhodata(:,8) + 1j* rhodata(:,11);
    rho13 = rhodata(:,9) + 1j* rhodata(:,15);

    rho20 = rhodata(:,4) - 1j* rhodata(:,10);
    rho21 = rhodata(:,8) - 1j* rhodata(:,11);
    rho22 = rhodata(:,12);
    rho23 = rhodata(:,13) + 1j* rhodata(:,16);

    rho30 = rhodata(:,5) - 1j* rhodata(:,14);
    rho31 = rhodata(:,9) - 1j* rhodata(:,15);
    rho32 = rhodata(:,13) - 1j* rhodata(:,16);
    rho33 = rhodata(:,17);

    % empty container for outputs
    fidelity = zeros(length(time), 1);  % fidelity of the state
    phi1_opt = zeros(length(time), 1);
    phi2_opt = zeros(length(time), 1);

    % container for [phi1, phi2]
    phi_last = [pi/2,pi/2];
    
    % for all timepoints in the series
    for i = 1:length(time)
        
        % find the trance of the density matrix
        Tr_rho = rho00(i) + rho11(i) + rho22(i) + rho33(i);
        
        % put the values in matrix form and normalize it with the trace
        rho = ...
            [rho00(i), rho01(i), rho02(i), rho03(i);...
             rho10(i), rho11(i), rho12(i), rho13(i);...
             rho20(i), rho21(i), rho22(i), rho23(i);...
             rho30(i), rho31(i), rho32(i), rho33(i)]...
             /Tr_rho;
        
        % minimize the infideity, oneminusF with FMINSEARCH
        oneminusF = @(phi) 1 - calculate_fidelity(rho, phi(1), phi(2));
        [phi_opt, oneminusF_min, exitflag] = fminsearch(oneminusF, [phi_last(1), phi_last(2)]); 
        
        % check for convergence
        if exitflag ~= 1
            error_msg = sprintf('fminsearch did not converge')
        end
        
        % save the optimal values
        phi_last = phi_opt;
        phi1_opt(i) = phi_opt(1);
        phi2_opt(i) = phi_opt(2);
        fidelity(i) = 1 - oneminusF_min;
    end
    
    % print the time series of fidelity to file, "fidelity_####.txt"
    dlmwrite(['./results/fidelity_', regstr, '.txt'], [time,fidelity], ...
        'delimiter', '\t',...
        'precision', 10);

end


