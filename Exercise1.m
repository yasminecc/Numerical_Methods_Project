% Exercise 1 

% Given parameter d = 1

d = 1;

% Parameters
R = 5 + 0.2 * d;  % Radius of the tank
V = 40;           % volume (m^3)

% Function h_func
h_func = @(h) ((pi*(h^2))/3)*(3*R - h) - V;

% Convergence tolerance
e = 1e-6;

%% (a) Bisection Method

fprintf('Bisection Method:\n');

a0 = 0.1;   % Initial lower bound
b0 = 10.0;  % Initial upper bound
error_bisec = abs(b0 - a0);  % Initial error estimate

% Bisection loop
while error_bisec > e
    h_mid = (a0 + b0) / 2;
    if h_func(a0)*h_func(h_mid) < 0
        b0 = h_mid;
    else
        a0 = h_mid;
    end
    
    error_bisec = abs(b0 - a0);   % Update error estimate
end

h_bisec = (a0 + b0) / 2;

% Display results
fprintf('Estimated depth (h) = %.3f m\n\n', h_bisec);


%% (b) Secant Method

fprintf('Secant Method:\n');

% Initial conditions
h_1 = 0.5; % h_-1
h_0 = 1.0; % h_0
e = 1e-6;

error_secant = Inf;

while error_secant > e
    f_1 = h_func(h_1);
    f_0   = h_func(h_0);
    h_next = h_0 - f_0 * ((h_0 - h_1) / (f_0 - f_1));
    error_secant = abs(h_next - h_0);

    % Update for next iteration
    h_1 = h_0;
    h_0 = h_next;
end

h_secant = h_0;

fprintf('Estimated depth (h) = %.3f m\n\n', h_secant);

%% (b) Newton-Raphson Method

fprintf('Newton-Raphson Method:\n');

% Initial guess
h_newton = 1.0;
e = 1e-6;
error_newton = Inf;

% Derivative of h_func(h)
h_prime = @(h) pi*(2*R*h - h^2);

while error_newton > e
    h_next = h_newton - h_func(h_newton)/h_prime(h_newton);
    error_newton = abs(h_next - h_newton);
    h_newton = h_next;
end

fprintf('Estimated depth (h) = %.3f m\n\n', h_newton);

