% Define the evaluation function for Set A (circle)
eval_A = @(x) norm(x) - 1;  % Distance of point x from the circle

% Define the projection function for Set A (circle)
prox_A = @(x,~) x / norm(x);  % Project x onto the circle's boundary

% Define the evaluation function for Set B (square)
eval_B = @(x) max(abs(x)) - 1;  % Distance of point x from the square

% Define the projection function for Set B (square)
prox_B = @(x,~) min(max(x, -1), 1);  % Project x onto the square's boundary

% Create the function structures for each set
F_A.eval = eval_A;
F_A.prox = prox_A;

F_B.eval = eval_B;
F_B.prox = prox_B;

% Combine the function structures into an array
F = {F_A, F_B};

% Set the starting point
x_0 = [0.5; 0.5];  % Starting point inside the circle but outside the square

% Call the pocs function to find the solution
[sol, info] = pocs(x_0, F);

% Plotting
figure;
hold on;
grid on;
axis equal;

% Plot Set A (circle)
theta = linspace(0, 2*pi, 100);
circle_x = cos(theta);
circle_y = sin(theta);
plot(circle_x, circle_y, 'b-', 'LineWidth', 2);
title('POCS Algorithm - Projection onto Convex Sets');
xlabel('x');
ylabel('y');

% Plot Set B (square)
square_x = [-1 -1 1 1 -1];
square_y = [-1 1 1 -1 -1];
plot(square_x, square_y, 'r-', 'LineWidth', 2);

% Plot the starting point
plot(x_0(1), x_0(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
text(x_0(1), x_0(2), 'Start', 'VerticalAlignment', 'bottom');

% Plot the solution
plot(sol(1), sol(2), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
text(sol(1), sol(2), 'Solution', 'VerticalAlignment', 'top');

% Display convergence information
text(-1, 1.5, ['Iterations: ' num2str(info.iter)], 'Color', 'b');
text(-1, 1.3, ['Final Norm: ' num2str(info.final_eval)], 'Color', 'b');
text(-1, 1.1, ['Stopping Criterion: ' info.crit], 'Color', 'b');

% Set axis limits
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);

% Legend
legend('Set A (circle)', 'Set B (square)', 'Starting Point', 'Solution');
