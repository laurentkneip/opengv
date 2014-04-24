%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% central case -> only one camera
cam_number = 4;
% let's only get 6 points, and generate new ones in each iteration
pt_number = 50;
% noise test, so no outliers
outlier_fraction = 0.0;
% repeat 1000 iterations
iterations = 1000;

% The algorithms we want to test
algorithms = { 'gp3p'; 'gpnp' };
% This defines the number of points used for every algorithm
indices = { [1, 2, 3]; [1:1:50] };
% The name of the algorithms on the plots
names = { 'GP3P'; 'GPnP (50pts)' };

% The noise in this experiment
noise = 1.0;

%% Run the benchmark

%prepare the overall result arrays
num_algorithms = size(algorithms,1);
execution_times = zeros(num_algorithms,iterations);
counter = 0;

for i=1:iterations

    % generate experiment
    [points,v,t,R] = create2D3DExperiment(pt_number,cam_number,noise,outlier_fraction);
    [t_perturbed,R_perturbed] = perturb(t,R,0.01);
    T_perturbed = [R_perturbed,t_perturbed];

    for a=1:num_algorithms
        tic
        T = opengv_donotuse(algorithms{a},indices{a},points,v,T_perturbed);
        execution_times(a,i) = toc/20.0;
    end

    counter = counter + 1;
    if counter == 100
        counter = 0;
        display(['Iteration ' num2str(i) ' of ' num2str(iterations) '(noise level ' num2str(noise) ')']);
    end
end

%% Plot the results

hist(execution_times')
legend(names,'Location','NorthWest')
xlabel('execution times [s]')
grid on

mean(execution_times')
