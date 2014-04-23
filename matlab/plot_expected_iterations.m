%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% The algorithms we want to test
algorithms = [ 6; 8; 17 ];
% The name of the algorithms in the final plots
names = { '6pt'; 'ge (8pt)'; '17pt'};

% The main experiment parameters
min_outlier_fraction = 0.05;%0.05;
max_outlier_fraction = 0.25;
outlier_fraction_step = 0.025;

p = 0.99;

%% Run the benchmark

%prepare the overall result arrays
number_outlier_fraction_levels = round((max_outlier_fraction - min_outlier_fraction) / outlier_fraction_step + 1);
num_algorithms = size(algorithms,1);
expected_number_iterations = zeros(num_algorithms,number_outlier_fraction_levels);
outlier_fraction_levels = zeros(1,number_outlier_fraction_levels);

%Run the experiment
for n=1:number_outlier_fraction_levels

    outlier_fraction = (n - 1) * outlier_fraction_step + min_outlier_fraction;
    outlier_fraction_levels(1,n) = outlier_fraction;
    display(['Analyzing outlier fraction level: ' num2str(outlier_fraction)])

    %Now compute the mean and median value of the error for each algorithm
    for a=1:num_algorithms
        expected_number_iterations(a,n) = log(1-p)/log(1-(1-outlier_fraction)^(algorithms(a,1)));
    end
    
end

%% Plot the results

figure(1)
plot(outlier_fraction_levels,expected_number_iterations,'LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('outlier fraction')
ylabel('expected number iterations')
grid on