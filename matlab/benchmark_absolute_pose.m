%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% central case -> only one camera
cam_number = 1;
% let's only get 6 points, and generate new ones in each iteration
pt_number = 6;
% noise test, so no outliers
outlier_fraction = 0.0;
% repeat 5000 iterations per noise level
iterations = 5000;

% The algorithms we want to test
algorithms = { 'p3p_kneip'; 'p3p_gao'; 'epnp'; 'abs_nonlin_central'; 'upnp'; 'upnp' };
% This defines the number of points used for every algorithm
indices = { [1, 2, 3]; [1, 2, 3]; [1, 2, 3, 4, 5, 6]; [1, 2, 3, 4, 5, 6]; [1, 2, 3, 4, 5, 6]; [1, 2, 3] };
% The name of the algorithms on the plots
names = { 'P3P (Kneip)'; 'P3P (Gao)'; 'EPnP'; 'nonlinear optimization'; 'UPnP'; 'UPnP (minimal)' };

% The maximum noise to analyze
max_noise = 5.0;
% The step in between different noise levels
noise_step = 0.1;

%% Run the benchmark

%prepare the overall result arrays
number_noise_levels = max_noise / noise_step + 1;
num_algorithms = size(algorithms,1);
mean_position_errors = zeros(num_algorithms,number_noise_levels);
mean_rotation_errors = zeros(num_algorithms,number_noise_levels);
median_position_errors = zeros(num_algorithms,number_noise_levels);
median_rotation_errors = zeros(num_algorithms,number_noise_levels);
noise_levels = zeros(1,number_noise_levels);

%Run the experiment
for n=1:number_noise_levels

    noise = (n - 1) * noise_step;
    noise_levels(1,n) = noise;
    display(['Analyzing noise level: ' num2str(noise)])
    
    position_errors = zeros(num_algorithms,iterations);
    rotation_errors = zeros(num_algorithms,iterations);
    
    counter = 0;
    
    for i=1:iterations
        
        % generate experiment
        [points,v,t,R] = create2D3DExperiment(pt_number,cam_number,noise,outlier_fraction);
        [t_perturbed,R_perturbed] = perturb(t,R,0.01);
        T_perturbed = [R_perturbed,t_perturbed];
        T_gt = [R,t];
        
        % run all algorithms
        for a=1:num_algorithms            
            T = opengv(algorithms{a},indices{a},points,v,T_perturbed);
            [position_error, rotation_error] = evaluateTransformationError( T_gt, T );
            position_errors(a,i) = position_error;
            rotation_errors(a,i) = rotation_error;
        end
        
        counter = counter + 1;
        if counter == 100
            counter = 0;
            display(['Iteration ' num2str(i) ' of ' num2str(iterations) '(noise level ' num2str(noise) ')']);
        end        
    end

    %Now compute the mean and median value of the error for each algorithm
    for a=1:num_algorithms
        mean_position_errors(a,n) = mean(position_errors(a,:));
        median_position_errors(a,n) = median(position_errors(a,:));
        mean_rotation_errors(a,n) = mean(rotation_errors(a,:));
        median_rotation_errors(a,n) = median(rotation_errors(a,:));
    end
    
end

%% Plot the results

figure(1)
plot(noise_levels',mean_rotation_errors','LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('noise level [pix]')
ylabel('mean rot. error [rad]')
grid on

figure(2)
plot(noise_levels',median_rotation_errors','LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('noise level [pix]')
ylabel('median rot. error [rad]')
grid on

figure(3)
plot(noise_levels',mean_position_errors','LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('noise level [pix]')
ylabel('mean pos. error [m]')
grid on

figure(4)
plot(noise_levels',median_position_errors','LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('noise level [pix]')
ylabel('median pos. error [m]')
grid on