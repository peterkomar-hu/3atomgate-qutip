function Fsq = calculate_fidelity(rho, phi1, phi2)
% this function finds the square overlap between rho and 
% 0.5 * ( |00> + e^(i*phi2)*|01> + e^(i*phi1)*|10> - e^(i*(phi1 + phi2))|11> )
    
    % target state
    psi = 0.5* [...
        1; ...
        exp(1j * phi1); ...
        exp(1j * phi2); ...
        -exp(1j* (phi1 + phi2))];
    
    % square overlap
    Fsq = real(psi' * rho * psi);

end