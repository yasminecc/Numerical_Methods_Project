%% Exercise 2

d = 1;

% Parameters
m0 = 1; 
I = 1;    
R = 1;   

% Coordinates of Point R
x = 0.5;
y = 0.2 + 0.3*d;
z = 0.5;

% Constant factor in the integral
const_factor = (m0 * I) / (4 * pi);

% Limits of integration
a = 0;
b = 2*pi;

integrand = @(t) [ ...
    (R*z*cos(t)), ...
    (R*z*sin(t)), ...
    (R^2 - R*(x*cos(t) + y*sin(t))) ] ./ ...
    ( ((x - R*cos(t)).^2 + (y - R*sin(t)).^2 + z^2).^(3/2) );

%% (a) Trapezoidal Rule 

n_trap = 15;
t_trap = linspace(a, b, n_trap+1);
f_trap = zeros(n_trap+1, 3);
for i = 1:n_trap+1
    f_trap(i,:) = integrand(t_trap(i));
end

h_trap = (b - a) / n_trap;
Bx_trap = (h_trap/2) * ( f_trap(1,1) + 2*sum(f_trap(2:end-1,1)) + f_trap(end,1) );
By_trap = (h_trap/2) * ( f_trap(1,2) + 2*sum(f_trap(2:end-1,2)) + f_trap(end,2) );
Bz_trap = (h_trap/2) * ( f_trap(1,3) + 2*sum(f_trap(2:end-1,3)) + f_trap(end,3) );

B_trap = const_factor * [Bx_trap, By_trap, Bz_trap];
B_trap_mag = sqrt(sum(B_trap.^2));

fprintf('Trapezoidal: |B| = %.3e\n\n', B_trap_mag);

%% (b) Simpson’s 1/3 Rule 

n_simp = 9;
t_simp = linspace(a, b, n_simp);
f_simp = zeros(n_simp,3);
for i = 1:n_simp
    f_simp(i,:) = integrand(t_simp(i));
end

h_simp = (b - a) / (n_simp - 1);

% Indices for Simpson's rule:

odd_ind = 2:2:n_simp-1;   % 2,4,6,8
even_ind = 3:2:n_simp-2; % 3,5,7

Bx_simp = (h_simp/3) * ( f_simp(1,1) + f_simp(end,1) + 4*sum(f_simp(odd_ind,1)) + 2*sum(f_simp(even_ind,1)) );
By_simp = (h_simp/3) * ( f_simp(1,2) + f_simp(end,2) + 4*sum(f_simp(odd_ind,2)) + 2*sum(f_simp(even_ind,2)) );
Bz_simp = (h_simp/3) * ( f_simp(1,3) + f_simp(end,3) + 4*sum(f_simp(odd_ind,3)) + 2*sum(f_simp(even_ind,3)) );

B_simp = const_factor * [Bx_simp, By_simp, Bz_simp];
B_simp_mag = sqrt(sum(B_simp.^2));

fprintf('Simpson''s: |B| = %.3e\n\n', B_simp_mag);

%% (c) Gauss-Legendre Quadrature 

% Gauss-Legendre points and weights for n=5
gauss_points = [-0.9061798459, -0.5384693101, 0, 0.5384693101, 0.9061798459];
gauss_weights = [0.2369268851, 0.4786286705, 0.5688888889, 0.4786286705, 0.2369268851];

% Transforming to interval [0, 2π]
t_points = 0.5 * (gauss_points + 1) * 2 * pi;
weights = gauss_weights * pi;

% Function for the Biot-Savart integral
biot_savart_integrand = @(t) sqrt( ...
    (R * z .* cos(t)).^2 + (R * z .* sin(t)).^2 + ...
    (R^2 - R .* (x .* cos(t) + y .* sin(t))).^2 ...
) ./ ((x - R .* cos(t)).^2 + (y - R .* sin(t)).^2 + z^2).^(3/2);

% Compute the integral 
integrand_values = biot_savart_integrand(t_points);
B_gauss_mag = (m0 * I / (4 * pi)) * sum(weights .* integrand_values);

fprintf('Gauss-Legendre: |B| = %.3e\n\n', B_gauss_mag);




