% Define parameters
a = 1.0;
b = 0.5;
L = 50;  % Spatial region [0, L]
time_steps = 4;
num_points = 1000;  % Number of points in the spatial domain

% Define spatial domain
x = linspace(0, L, num_points);

% Initialize plot
figure;
hold on;

% Time evolution
for t = linspace(0, L, time_steps)
    psi = a * sqrt(2) * exp(1i * ((b * x * sqrt(2))/ 2 - (b^2 / 4 - a^2) * t)) .* sech(a * (x * sqrt(2) - b * t));
    
    % Plot amplitude (absolute value of psi)
    plot(x, abs(psi), 'DisplayName', ['t = ' num2str(t)]);
end

hold off;
xlabel('Spatial Domain');
ylabel('|ψNLS(x, t)|');
title('Time Evolution of ψNLS(x, t) between 0 and L');
legend('show');
grid on;
