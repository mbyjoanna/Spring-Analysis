% Vertical Case
% Define constants
m = 0.240; % Mass of phone in kg
num_trials = 6;
spring_constants = zeros(1, num_trials);
damping_coefficients = zeros(1, num_trials);
decay_rates = zeros(1, num_trials);
damping_ratios = zeros(1, num_trials);
natural_frequencies = zeros(1, num_trials);

figure; hold on; % Create a figure to plot all trials together
xlabel('Time (s)'); ylabel('Amplitude (m/s^2)');
title('Amplitude vs Time for All Trials');

for i = 1:num_trials
    % Load data from CSV file
    filename = sprintf('springv%d.csv', i); % Ensure file names match
    data = readmatrix(filename);  % Use readmatrix() instead of load()
    
    % Extract time and acceleration in y direction
    time = data(:, 1);
    acceleration = data(:, 3);
    
    % Plot amplitude vs time
    plot(time, acceleration, 'DisplayName', sprintf('Trial %d', i));
    
    % Find peaks using prominence filtering
    [pks, locs] = findpeaks(acceleration, time, 'MinPeakProminence', 0.1);
    
    if length(pks) < 2
        warning('Not enough significant peaks detected in trial %d', i);
        continue;
    end

    % Compute periods from peak locations
    periods = diff(locs);
    avg_period = mean(periods);
    omega_d = 2 * pi / avg_period;
    
    % Logarithmic decrement method for decay rate
    delta = log(pks(1:end-1) ./ pks(2:end));
    decay_rate = abs(mean(delta) / mean(periods)); % Ensure decay rate is positive
    
    % Calculate damping ratio
    zeta = decay_rate / sqrt(decay_rate^2 + omega_d^2);
    
    % Calculate natural frequency
    omega_n = omega_d / sqrt(1 - zeta^2);
    
    % Calculate spring constant
    k = omega_n^2 * m;
    
    % Calculate damping coefficient
    b = 2 * m * omega_n * zeta;
    
    % Store results
    spring_constants(i) = k;
    damping_coefficients(i) = b;
    decay_rates(i) = decay_rate;
    damping_ratios(i) = zeta;
    natural_frequencies(i) = omega_n;
end

legend; % Show legend for different trials

% Compute average values
avg_k = mean(spring_constants);
std_k = std(spring_constants); % Standard deviation of k
avg_b = mean(damping_coefficients);
avg_decay_rate = mean(decay_rates);
avg_zeta = mean(damping_ratios);
avg_omega_n = mean(natural_frequencies);

% Display results
fprintf('Average results for vertical orientation:\n');
fprintf('Average Spring Constant: %.4f N/m\n', avg_k);
fprintf('Standard Deviation of Spring Constant: %.4f N/m\n', std_k);
fprintf('Average Damping Coefficient: %.10f Ns/m\n', avg_b);
fprintf('Average Decay Rate: %.10f s⁻¹\n', avg_decay_rate);
fprintf('Average Damping Ratio: %.10f\n', avg_zeta);
fprintf('Average Natural Frequency: %.10f rad/s\n', avg_omega_n)
