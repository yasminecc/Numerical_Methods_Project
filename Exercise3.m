%% Exercise 3

% Given parameters
d = 1; 
kg_max = 0.3 + 0.1 * d;  
Y = 0.5;                 
Ks = 120;                
kd = 0.01;               
kr = 0.01;              
Tw = 20;             
S_in = 1000;              

% Initial guesses for X and S
X = 100; 
S = 0;   

% Convergence parameters
iterations = 5;

for i = 1:iterations
    % Compute function values and Jacobian
    [F, J] = sysEqu(X, S, kg_max, Y, Ks, kd, kr, Tw, S_in);
    
    if i == 1
        % (a) The determinant of on the first iteration
        detJ = det(J);
        fprintf('Determinant of Jacobian matrix on the first iteration = %.3f\n', detJ);
    end
    
    delta = -J \ F;
    
    % Update X and S
    X = X + delta(1);
    S = S + delta(2);
end

%(b) Final result after 5 iterations
fprintf('Final Result (X,S) after 5 iterations: X = %.3f gC/m^3, S = %.3f gC/m^3\n', X, S);

%% Function to compute the system of equations and Jacobian
function [F, J] = sysEqu(X, S, kg_max, Y, Ks, kd, kr, tau_w, Sin)
   
    % Computing function values
    f1 = (kg_max * Y * S / (Ks + S)) * X - (kd + kr + 1/tau_w) * X;
    f2 = -(kg_max * S / (Ks + S)) * X + kd * X + (1/tau_w) * (Sin - S);
    F = [f1; f2];
    
    % Computing partial derivatives 
    dfdX1 = (kg_max * Y * S / (Ks + S)) - (kd + kr + 1/tau_w);
    dfdS1 = (kg_max * Y * X * Ks) / ((Ks + S)^2);
    dfdX2 = -(kg_max * S / (Ks + S)) + kd;
    dfdS2 = -(kg_max * X * Ks) / ((Ks + S)^2) - (1/tau_w);
    J = [dfdX1, dfdS1; dfdX2, dfdS2];
end
