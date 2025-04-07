% Exercise 6

% Parameters
z0 = 100; 
v0 = 60; 
g = 9.81; 
m = 90; 
d = 1; 
c = 5 + 0.8 * d;

% Newton-Raphson Parameters
e = 1e-6; % Convergence tolerance
i_max = 5; % Maximum number of iterations
t_guess = 0; % Initial guess for peak time

% Define Velocity and its Derivative
velocity = @(t) (v0 + (m * g) / c) * exp(-c * t / m) - (m * g) / c;
velocity_derivative = @(t) -(c / m) * (v0 + (m * g) / c) * exp(-c * t / m);

% Newton-Raphson Method to Find Peak Time
t_peak = t_guess;
for iter = 1:i_max
    t_new = t_peak - velocity(t_peak) / velocity_derivative(t_peak); % Update time
    if abs(t_new - t_peak) < e
        break; 
    end
    t_peak = t_new; % Update the estimate
end

% Compute Peak Altitude
z_peak = z0 + (m / c) * (v0 + (m * g) / c) * (1 - exp(-c * t_peak / m)) - (m * g / c) * t_peak;

% Display Results
fprintf('Time at peak altitude: %.3f seconds\n', t_peak);
fprintf('Peak altitude: %.3f meters\n', z_peak);

% Altitude as a Function of Time
t_vals = linspace(0, 2 * t_peak, 100); % Time range for plotting
altitude = @(t) z0 + (m / c) * (v0 + (m * g) / c) * (1 - exp(-c * t / m)) - (m * g / c) * t;
z_vals = altitude(t_vals);

% Plot Altitude vs Time
figure;
plot(t_vals, z_vals, 'm-', 'LineWidth', 1, 'DisplayName', 'Altitude');
hold on;
plot(t_peak, z_peak, 'b*', 'MarkerSize', 15, 'DisplayName', 'Peak Altitude');
xlabel('Time (s)');
ylabel('Altitude (m)');
title('Altitude vs Time');
grid on;
legend show;
