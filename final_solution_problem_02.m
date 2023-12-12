
num_collisions = 20; % Number of collisions to simulate
N = 6000;         % Number of spatial nodes
L = 800;          % Length of spatial domain
dx = L / N;       % Spatial step size
x = 0:dx:L;       % Spatial grid
% Initial parameters for the first soliton
a = 0.1;          % Amplitude parameter
b = 3;            % Frequency parameter
x0 = 100;         % Initial position of the first soliton
gif_frames = cell(1, num_collisions);
% Initial parameters for the second soliton
a1 = 0.2 * a;     % Amplitude parameter for the second soliton
b1 = -3;          % Frequency parameter for the second soliton
x1 = 700;         % Initial position of the second soliton

for i = 1:num_collisions
    % Calculate the wave function for each iteration
    psi = a * sqrt(2) * exp(1i * b * x / 2) .* sech(a * (x - x0)) + a1 * sqrt(2) * exp(1i * b1 * x / 2) .* sech(a1 * (x - x1));
    
    % For demonstration, let's assume a shift towards each other in each iteration
    shift_amount = 5 * i; % Increase the shift amount towards each other
    
    % Shift the positions towards each other
    x0 = x0 + shift_amount / 2; % Shift the first soliton leftwards
    x1 = x1 - shift_amount / 2; % Shift the second soliton rightwards
    
    % Plot the absolute value of psi after each iteration
    psi_abs = abs(psi); % Take the absolute value of psi
    
    % Create masks based on amplitude for two solitons
    mask1 = psi_abs >= a * sqrt(2); % Mask for soliton 1
    mask2 = psi_abs < a * sqrt(2);  % Mask for soliton 2
    
    % Plot the absolute value of psi after each iteration
    
    figure(); % Create an invisible figure
    plot(x(mask1), psi_abs(mask1), 'b'); % Plot soliton 1 with blue color
    hold on;
    plot(x(mask2), psi_abs(mask2), 'r'); % Plot soliton 2 with red color
    hold on;
    
    xlabel('Position');
    ylabel('|\psi(x,t)|');
    title(['Soliton-Soliton Collison after Iteration ', num2str(i)]);
    
    % Save the plot as an image within the loop
    filename = ['soliton_collision_' num2str(i) '.png']; % Define the filename
    saveas(gcf, filename); % Save the current figure
    
    % Capture the frame for the GIF
    gif_frames{i} = getframe(gcf);
    close(gcf); % Close the figure after capturing the frame
end

% % % Create a GIF from captured frames %this part need to be uncommented to
% % % save the GIF figure

% filename = 'soliton_collision.gif';
% for i = 1:num_collisions
%     frame = gif_frames{i}.cdata;
%     [imind, cm] = rgb2ind(frame, 256);
%     if i == 1
%         imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.5);
%     else
%         imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.5);
%     end
% end