function K = deviation(C, Omega, kappa, gamma, Delta_E, Delta_e)
% this function calculates the total deviation of the numerically
% determined decay rates of the three sectors: 0, 1, or 2 qubits in "1"

% its output is a single number, the logarithm of sum of square deviations
% the inputs are
% C: the cooperativity
% kappa: cavity decay rate
% gamma: excited atomic state decay rate
% Delta_E = Delta1 : detuning of the control atom's excited state
% Delta_e = Delta2 : detuning of the qubit atoms' excited state

% detmine the atom-cavity coupling from C
g = sqrt(C * kappa * gamma);

% Hamiltonians of the three sectors:

% qubits in 00 (total = 0)
iH0 = [ 0 , 1i*Omega/2, 0 ; ...
        1i*Omega/2, 1i*Delta_E + gamma/2, 1i*g ; ...
        0, 1i*g , kappa/2];
    
% qubits in 01 or 10 (total = 1)
iH1 = [ 0 , 1i*Omega/2, 0 , 0; ...
        1i*Omega/2, 1i*Delta_E + gamma/2, 1i*g , 0; ...
        0, 1i*g , kappa/2, 1i*g; ...
        0 , 0 , 1i*g , 1i*Delta_e + gamma/2];
    
% qubits in 11 (total = 2)
iH2 = [ 0 , 1i*Omega/2, 0 , 0, 0; ...
        1i*Omega/2, 1i*Delta_E + gamma/2, 1i*g , 0, 0; ...
        0, 1i*g , kappa/2, 1i*g, 1i*g; ...
        0 , 0 , 1i*g , 1i*Delta_e + gamma/2, 0; ...
        0 , 0 , 1i*g , 0 , 1i*Delta_e + gamma/2];
    
% calculate eigenvalues
lambda0 = eig(iH0);
lambda1 = eig(iH1);
lambda2 = eig(iH2);

% determine decay rates
kappa0 = 2 * min(real(lambda0));
kappa1 = 2 * min(real(lambda1));
kappa2 = 2 * min(real(lambda2));

% return the log of square sum of pariwise differences, as deviation
K = log( (kappa0 - kappa1)^2 + (kappa0 - kappa2)^2 + (kappa1 - kappa2)^2 );

end




