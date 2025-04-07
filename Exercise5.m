% Exercise 5

% Parameters
G = 1.0;
d = 1;   
m = 0.1 * (1 + d); 
L = 1;   
U = 1;  

n = 10;        % Number of segments
h = L / n;     % Step size
y = linspace(0, L, n+1); % Grid points

% coefficient matrix
A = zeros(n-1, n-1);
b = -G / m * ones(n-1, 1);

% Fill the tridiagonal matrix
for i = 1:n-1
    if i > 1
        A(i, i-1) = 1 / h^2;
    end
    A(i, i) = -2 / h^2;
    if i < n-1
        A(i, i+1) = 1 / h^2;
    end
end

b(end) = b(end) - U / h^2; 

% Solve the linear system
u_fd = A \ b;
u_fd = [0; u_fd; U]; % boundary values at y=0 and y=L

% Exact solution
u_exact = (G / (2 * m)) * y .* (L - y) + (U / L) * y;

% Plot the results
figure;
plot(y, u_fd, 'o-', 'LineWidth', 1.25, 'DisplayName', 'Finite Difference');
hold on;
plot(y, u_exact, '-', 'LineWidth', 1.25, 'DisplayName', 'Exact Solution');
xlabel('y');
ylabel('u(y)');
title('Velocity Profile for Couette Flow');
legend('Location', 'best');
grid on;

% Print results in requested format
fprintf('grid point yi | calculated velocity\n');
for i = 1:length(y)
    fprintf('%.4f          %.4f\n', y(i), u_fd(i));
end
