%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% noncentral case
cam_number = 4;
% set maximum and minimum number of points per cam
pt_number_per_cam =  50;
% set maximum and minimum number of outliers
min_outlier_fraction = 0.1;
max_outlier_fraction = 0.2;
% repeat 10000 iterations
iterations = 1000;

% The name of the algorithms in the final plots
names = { 'Homogeneous'; 'Vanilla' };

% The noise in this experiment
noise = 0.5;

%% Run the benchmark

%prepare the overall result arrays
ransac_iterations = zeros(2,iterations);

%Run the RANSAC with homogeneous sampling
counter = 0;

for i=1:iterations

    %generate random outlier fraction
    outlier_fraction = rand() * (max_outlier_fraction - min_outlier_fraction) + min_outlier_fraction;
    
    % generate experiment        
    [v1,v2,cam_offsets,t,R] = createMulti2D2DExperiment(pt_number_per_cam,cam_number,noise,outlier_fraction);
    Out = opengv_experimental1( v1{1,1}, v1{2,1}, v1{3,1}, v1{4,1}, v2{1,1}, v2{2,1}, v2{3,1}, v2{4,1}, cam_offsets, 2 );
    ransac_iterations(1,i) = Out(1,5);

    counter = counter + 1;
    if counter == 100
        counter = 0;
        display(['Homogeneous sampling: Iteration ' num2str(i) ' of ' num2str(iterations)]);
    end        
end

%Run the RANSAC with vanilla sampling
counter = 0;

for i=1:iterations

    %generate random outlier fraction
    outlier_fraction = rand() * (max_outlier_fraction - min_outlier_fraction) + min_outlier_fraction;
    
    % generate experiment        
    [v1,v2,t,R] = create2D2DExperiment(pt_number_per_cam*cam_number,cam_number,noise,outlier_fraction);
    Out = opengv_experimental2( v1, v2, 2 );
    ransac_iterations(2,i) = Out(1,5);

    counter = counter + 1;
    if counter == 100
        counter = 0;
        display(['Vanilla sampling: Iteration ' num2str(i) ' of ' num2str(iterations)]);
    end        
end

%% Plot the results

figure(1)
hist(ransac_iterations')
[Y,X] = hist(ransac_iterations')
legend(names,'Location','NorthEast')
xlabel('number of iterations')
grid on