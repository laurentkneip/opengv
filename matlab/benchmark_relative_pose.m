%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% central case -> only one camera
cam_number = 1;
% Getting 10 points, and testing all algorithms with the respective number of points
pt_number = 10;
% noise test, so no outliers
outlier_fraction = 0.0;
% repeat 5000 tests per noise level
iterations = 5000;

% The algorithms we want to test
algorithms = { 'fivept_stewenius'; 'fivept_nister'; 'fivept_kneip'; 'sevenpt'; 'eightpt'; 'eigensolver'; 'rel_nonlin_central' };
% Some parameter that tells us what the result means
returns = [ 1, 1, 0, 1, 1, 0, 2 ]; % 1means essential matrix(ces) needing decomposition, %0 means rotation matrix(ces), %2 means transformation matrix
% This defines the number of points used for every algorithm
indices = { [1, 2, 3, 4, 5]; [1, 2, 3, 4, 5]; [1, 2, 3, 4, 5]; [1, 2, 3, 4, 5, 6, 7]; [1, 2, 3, 4, 5, 6, 7, 8]; [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]; [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] };
% The name of the algorithms in the final plots
names = { '5pt (Stewenius)'; '5pt (Nister)'; '5pt (Kneip)'; '7pt'; '8pt'; 'eigensolver (10pts)'; 'nonlin. opt. (10pts)' };

% The maximum noise to analyze
max_noise = 5.0;
% The step in between different noise levels
noise_step = 0.1;

%% Run the benchmark

%prepare the overall result arrays
number_noise_levels = max_noise / noise_step + 1;
num_algorithms = size(algorithms,1);
mean_rotation_errors = zeros(num_algorithms,number_noise_levels);
median_rotation_errors = zeros(num_algorithms,number_noise_levels);
noise_levels = zeros(1,number_noise_levels);

%Run the experiment
for n=1:number_noise_levels

    noise = (n - 1) * noise_step;
    noise_levels(1,n) = noise;
    display(['Analyzing noise level: ' num2str(noise)])
    
    rotation_errors = zeros(num_algorithms,iterations);
    
    counter = 0;
    
    validIterations = 0;
    
    for i=1:iterations
        
        % generate experiment        
        [v1,v2,t,R] = create2D2DExperiment(pt_number,cam_number,noise,outlier_fraction);
        [t_perturbed,R_perturbed] = perturb(t,R,0.01);
        T_perturbed = [R_perturbed,t_perturbed];
        R_gt = R;
        
        % run all algorithms
        allValid = 1;
        
        for a=1:num_algorithms
            Out = opengv(algorithms{a},indices{a},v1,v2,T_perturbed);
            
            if ~isempty(Out)
            
                if returns(1,a) == 1
                    temp = transformEssentials(Out);
                    Out = temp;
                end
                if returns(1,a) == 2
                    temp = Out(:,1:3);
                    Out = temp;
                end
            
                rotation_errors(a,validIterations+1) = evaluateRotationError( R_gt, Out );
                
            else
                
                allValid = 0;
                break;
                
            end
        end
        
        if allValid == 1
            validIterations = validIterations +1;
        end
        
        counter = counter + 1;
        if counter == 100
            counter = 0;
            display(['Iteration ' num2str(i) ' of ' num2str(iterations) '(noise level ' num2str(noise) ')']);
        end        
    end

    %Now compute the mean and median value of the error for each algorithm
    for a=1:num_algorithms
        mean_rotation_errors(a,n) = mean(rotation_errors(a,1:validIterations));
        median_rotation_errors(a,n) = median(rotation_errors(a,1:validIterations));
    end
    
end

%% Plot the results

figure(1)
plot(noise_levels,mean_rotation_errors,'LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('noise level [pix]')
ylabel('mean rot. error [rad]')
grid on

figure(2)
plot(noise_levels,median_rotation_errors,'LineWidth',2)
legend(names,'Location','NorthWest')
xlabel('noise level [pix]')
ylabel('median rot. error [rad]')
grid on
