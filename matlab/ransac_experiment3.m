%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% The name of the algorithms in the final plots
names = { '6pt'; 'ge (8pt)'; '17pt'};

% The main experiment parameters
min_outlier_fraction = 0.05;
max_outlier_fraction = 0.45;
outlier_fraction_step = 0.05;

%% Run the benchmark

%prepare the overall result arrays
number_outlier_fraction_levels = round((max_outlier_fraction - min_outlier_fraction) / outlier_fraction_step + 1);
num_algorithms = size(names,1);
mean_number_iterations = zeros(num_algorithms,number_outlier_fraction_levels);
mean_execution_times = zeros(num_algorithms,number_outlier_fraction_levels);
outlier_fraction_levels = zeros(1,number_outlier_fraction_levels);

%Run the experiment
for n=1:number_outlier_fraction_levels

    outlier_fraction = (n - 1) * outlier_fraction_step + min_outlier_fraction;
    outlier_fraction_levels(1,n) = outlier_fraction;
    display(['Analyzing outlier fraction level: ' num2str(outlier_fraction)])
    
    clear number_iterations
    clear execution_times
    temp_file_name1 = ['number_iterations_' num2str(outlier_fraction) '.mat'];
    temp_file_name2 = ['execution_times_' num2str(outlier_fraction) '.mat'];
    load(temp_file_name1)
    load(temp_file_name2)
    
    %Now compute the mean and median value of the error for each algorithm
    for a=1:num_algorithms
        mean_number_iterations(a,n) = mean(number_iterations(a,:));
        mean_execution_times(a,n) = mean(execution_times(a,:));
    end
end

%% Plot the results
 
figure(1)
plot(outlier_fraction_levels,mean_number_iterations,'LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('outlier fraction')
ylabel('mean number iterations')
axis([0.05 0.25 0 1500])
grid on

figure(2)
hold on
plot(outlier_fraction_levels(1,1:6),mean_execution_times(1,1:6),'LineWidth',2)
plot(outlier_fraction_levels,mean_execution_times(2:3,:),'LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('outlier fraction')
ylabel('mean execution time [s]')
axis([0.05 0.45 0 40])
grid on