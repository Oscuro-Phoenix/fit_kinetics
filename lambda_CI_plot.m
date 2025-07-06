% Example data
lambda = [4.11, 5.07, 4.69, 4.65, 4.9]; % Values of lambda
lambda_error = [0.26, 0.6, 0.38, 0.21, 0.22]; % Confidence intervals (plus-minus values)
x_labels = {'LCO LiClO4', 'LCO LiPF6', 'NMC111 LiClO4', 'NMC111 LiPF6', 'NMC811'}; % X-axis labels (strings)

% Create figure
figure('Position', [100, 100, 600, 400]); % Adjust figure size for publication

% Plot lambda values with error bars
errorbar(lambda, lambda_error, 'o', 'MarkerSize', 10, 'LineWidth', 2, 'CapSize', 10, 'Color', [0, 0.5, 0]);

% Set the x-axis labels and other axis properties
set(gca, 'XTick', 1:length(lambda), 'XTickLabel', x_labels, 'FontSize', 12, 'LineWidth', 1.5);
xtickangle(45); % Rotate x-axis labels for better readability

% Set Y-axis limits to provide clear spacing around the data points
ylim([min(lambda - lambda_error) - 0.5, max(lambda + lambda_error) + 0.5]);
padding = 0.5; % Amount of padding to add (adjust as needed)
xlim([0.5 - padding, length(lambda) + 0.5 + padding]); % Add padding to both ends of x-axis

% Add axis labels with appropriate font size
xlabel('Electrode-Electrolyte Combination', 'FontSize', 14);
ylabel('$\tilde{\lambda}$', 'Interpreter', 'latex', 'FontSize', 14);

% Title for the plot
title('Dimensionless \lambda with Confidence Intervals', 'FontSize', 16);

% Add grid lines for clarity
grid on;

% Make sure that the plot looks clean by removing unnecessary borders
set(gca, 'box', 'on');
% Optional: Save the figure as a high-resolution image for publication
saveas(gcf, 'lambda_with_error_bars.png'); % Save as PNG (you can change format to EPS, PDF, etc.)

% Optional: Adjust the figure to be publication quality
set(gcf, 'Color', 'w'); % Set background color to white
