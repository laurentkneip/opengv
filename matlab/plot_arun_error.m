%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% noncentral case
cam_number = 4;
% Getting 17 points, and testing all algorithms with the respective number of points
pt_number = 8;
% noise test, so no outliers
outlier_fraction = 0.0;
% repeat 1000 tests per noise level
iterations = 10000;

% The algorithms we want to test
algorithms = { 'ge2' };
% This defines the number of points used for every algorithm
indices = { [1:1:8] };
% The name of the algorithms in the final plots
names = { 'arun (8pt)' };

% The maximum noise to analyze
noise = 0.5;

%% Run the benchmark

%Run the experiment
    
rotation_errors = zeros(1,iterations);
    
counter = 0;
    
for i=1:iterations

    % generate experiment        
    [v1,v2,t,R] = create2D2DOmniExperiment(pt_number,cam_number,noise,outlier_fraction);
    [t_perturbed,R_perturbed] = perturb(t,R,0.01);
    T_perturbed = [R_perturbed,t_perturbed];
    T_init = [eye(3),zeros(3,1)];
    T_gt = [R,t];
        
    Out = opengv(algorithms{1},indices{1},v1,v2,T_perturbed);

    rotation_error = evaluateRotationError( R, Out(1:3,1:3) );
    rotation_errors(1,i) = rotation_error;
	
        
    counter = counter + 1;
    if counter == 100
        counter = 0;
        display(['Iteration ' num2str(i) ' of ' num2str(iterations) '(noise level ' num2str(noise) ')']);
    end        
end

%% Plot the results
hist(rad2deg(rotation_errors))
xlabel('rotation error [deg]')
ylabel('occurence')
grid on