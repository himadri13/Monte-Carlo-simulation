clear all; clc;

% Simulation window parameters (annular disc)
rMin = 1; % Inner radius
rMax = 2; % Outer radius
areaTotal = pi * (rMax^2 - rMin^2); % Area of the annular region

% User inputs
alpha = input('Enter the path-loss exponent (alpha): ');
power = input('Enter the power of the points: ');
numSimulations = 500; % Reduce number of Monte Carlo simulations to 500 for more visible variation

% Mean densities to be tested (from 10 to 30 with step of 2)
meanDensities = 10:2:30;

% Array to store the mean interference for each mean density (simulated)
meanInterferenceValuesSimulated = zeros(length(meanDensities), 1);

% Array to store the theoretical mean interference for each mean density
meanInterferenceValuesTheoretical = zeros(length(meanDensities), 1);

% For plotting individual simulation runs (optional)
interferenceSamples = zeros(numSimulations, length(meanDensities));

% Loop over the different mean densities
for idx = 1:length(meanDensities)
    lambda = meanDensities(idx); % Mean density (intensity) for this iteration

    % Initialize array to store interference values for this density (simulated)
    interferenceValues = zeros(numSimulations, 1);

    % Monte Carlo simulation
    for sim = 1:numSimulations
        % Simulate Poisson point process
        numbPoints = poissrnd(areaTotal * lambda); % Poisson number of points

        % Generate random points in polar coordinates
        theta = 2 * pi * rand(numbPoints, 1); % Uniform random angles between 0 and 2*pi
        r = sqrt((rMax^2 - rMin^2) * rand(numbPoints, 1) + rMin^2); % Random radius

        % Compute interference from each point at the center (origin)
        interference = sum(power ./ (r .^ alpha)); % Sum of interferences

        % Store the total interference for this simulation
        interferenceValues(sim) = interference;
    end

    % Estimate mean interference for this mean density (simulated)
    meanInterferenceValuesSimulated(idx) = mean(interferenceValues);

    % Store individual samples to plot later if needed
    interferenceSamples(:, idx) = interferenceValues;

    % Compute the theoretical mean interference for this mean density
    if alpha ~= 2
        % Use the general formula when alpha != 2
        meanInterferenceValuesTheoretical(idx) = (2 * pi * lambda * power) * ...
            (((rMax)^(2 - alpha)) - (rMin)^(2 - alpha)) * (1 / (2 - alpha));
    else
        % Use the logarithmic formula when alpha == 2
        meanInterferenceValuesTheoretical(idx) = 2 * pi * lambda * power * log(rMax / rMin);
    end
end

% Plotting the results (simulated vs theoretical)
figure;
hold on;

% Plot simulated mean interference
plot(meanDensities, meanInterferenceValuesSimulated, '-o', 'LineWidth', 2, 'MarkerSize', 8, 'DisplayName', 'Simulated');

% Plot theoretical mean interference
plot(meanDensities, meanInterferenceValuesTheoretical, '-s', 'LineWidth', 2, 'MarkerSize', 8, 'Color', 'r', 'DisplayName', 'Theoretical');

% Plot some individual simulation runs to show variability (optional)
for i = 1:10
    plot(meanDensities, interferenceSamples(randi(numSimulations), :), ':', 'Color', [0.6 0.6 0.6], 'DisplayName', 'Individual Sim.');
end

hold off;

% Labeling the plot
xlabel('Mean Density (\lambda)', 'FontSize', 12);
ylabel('Mean Interference at the Origin', 'FontSize', 12);
title('Simulated vs Theoretical Mean Interference', 'FontSize', 14);
legend('show'); % Show legend
grid on;
