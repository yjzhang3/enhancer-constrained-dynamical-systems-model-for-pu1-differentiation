% 1. Define the system of differential equations
dxdt = @(x,y) y - x.^2 + 1;
dydt = @(x,y) y + x.^2 - 1;

figure; hold on;

% 2. Plot the x-nullcline (dx/dt = 0)
fimplicit(@(x,y) dxdt(x,y), [-3 3 -3 3], 'b', 'LineWidth', 2);

% 3. Plot the y-nullcline (dy/dt = 0)
fimplicit(@(x,y) dydt(x,y), [-3 3 -3 3], 'r', 'LineWidth', 2);

% 4. Add a quiver plot for the phase portrait/vector field
[X, Y] = meshgrid(-3:0.5:3, -3:0.5:3);
U = dxdt(X, Y);
V = dydt(X, Y);
quiver(X, Y, U, V, 'k');

% Formatting the plot
xlabel('x');
ylabel('y');
title('Phase Portrait and Nullclines');
legend('dx/dt = 0 (x-nullcline)', 'dy/dt = 0 (y-nullcline)', 'Vector Field');
axis([-3 3 -3 3]);
grid on;
