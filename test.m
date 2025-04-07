%% Exercise 4 

% Parameters
m = 1; 
a = 5; 
k = 6; 
d = 1; 
F0 = 2.1 + 0.8 * d; 
w = 0.5; 
h = 0.125; 
t_upper = 15; 

% Define the force function as a nested function for clarity
force = @(t, v, u) (F0 * sin(w * t) - a * abs(v) * v - k * u) / m;

% Time vector
t = 0:h:t_upper;
n = length(t);

% Initial conditions
initial_u = 1; % Initial displacement at t=0
initial_v = 0; % Initial velocity at t=0

%% Euler's Method
[x_euler, v_euler] = eulerMethod(force, t, h, initial_u, initial_v);

% Display Euler Results
disp('Euler Results:');
fprintf('Displacement at t = %.1f: %.3f m\n', t_upper, x_euler(end));

%% Runge-Kutta (RK4) Method
[x_rk, v_rk] = rungeKuttaMethod(force, t, h, initial_u, initial_v);

% Display RK4 Results
disp('Runge-Kutta Results:');
fprintf('Displacement at t = %.1f: %.3f m\n', t_upper, x_rk(end));

%% Plot Results
figure;

% Euler's Method subplot
subplot(2, 1, 1);
plot(t, x_euler, 'b-', 'LineWidth', 1.25, 'DisplayName', 'Displacement');
hold on;
plot(t, v_euler, 'm--', 'LineWidth', 1.25, 'DisplayName', 'Velocity');
xlabel('Time (sec)');
ylabel('Displacement (m) / Velocity (m/s)');
title('Euler''s Method');
legend('Location', 'best');
grid on;

% Runge-Kutta Method subplot
subplot(2, 1, 2);
plot(t, x_rk, 'b-', 'LineWidth', 1.25, 'DisplayName', 'Displacement');
hold on;
plot(t, v_rk, 'm--', 'LineWidth', 1.25, 'DisplayName', 'Velocity');
xlabel('Time (sec)');
ylabel('Displacement (m) / Velocity (m/s)');
title('Fourth-order Runge-Kutta Method');
legend('Location', 'best');
grid on;

%% Functions

function [x_vals, v_vals] = eulerMethod(forceFunc, t, h, u0, v0)
    % Preallocate arrays
    x_vals = zeros(size(t));
    v_vals = zeros(size(t));
    
    % Initial conditions
    x_vals(1) = u0;
    v_vals(1) = v0;
    
    for i = 1:length(t)-1
        du_dt = v_vals(i);
        dv_dt = forceFunc(t(i), v_vals(i), x_vals(i));
        
        x_vals(i+1) = x_vals(i) + h * du_dt;
        v_vals(i+1) = v_vals(i) + h * dv_dt;
    end
end

function [x_vals, v_vals] = rungeKuttaMethod(forceFunc, t, h, u0, v0)
    % Preallocate arrays
    x_vals = zeros(size(t));
    v_vals = zeros(size(t));
    
    % Initial conditions
    x_vals(1) = u0;
    v_vals(1) = v0;
    
    for i = 1:length(t)-1
        % k1 terms
        k1u = v_vals(i);
        k1v = forceFunc(t(i), v_vals(i), x_vals(i));
        
        % k2 terms
        k2u = v_vals(i) + 0.5 * h * k1v;
        k2v = forceFunc(t(i) + 0.5 * h, v_vals(i) + 0.5*h*k1v, x_vals(i) + 0.5*h*k1u);
        
        % k3 terms
        k3u = v_vals(i) + 0.5 * h * k2v;
        k3v = forceFunc(t(i) + 0.5 * h, v_vals(i) + 0.5*h*k2v, x_vals(i) + 0.5*h*k2u);
        
        % k4 terms
        k4u = v_vals(i) + h * k3v;
        k4v = forceFunc(t(i) + h, v_vals(i) + h*k3v, x_vals(i) + h*k3u);
        
        x_vals(i+1) = x_vals(i) + (h/6)*(k1u + 2*k2u + 2*k3u + k4u);
        v_vals(i+1) = v_vals(i) + (h/6)*(k1v + 2*k2v + 2*k3v + k4v);
    end
end
