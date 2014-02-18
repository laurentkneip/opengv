%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% central case -> only one camera
cam_number = 4;
% Getting 17 points, and testing all algorithms with the respective number of points
pt_number = 17;
% noise test, so no outliers
outlier_fraction = 0.0;
% repeat 1000 iterations
iterations = 1000;

% The algorithms we want to test
algorithms = { 'sixpt'; 'ge'; 'seventeenpt' };
% This defines the number of points used for every algorithm
indices = { [1:1:6]; [1:1:8]; [1:1:17] };
% The name of the algorithms in the final plots
names = { '6pt';'ge (8pt)'; '17pt' };

% The noise in this experiment
noise = 0.5;

%% Run the benchmark

%prepare the overall result arrays
num_algorithms = size(algorithms,1);
execution_times = zeros(num_algorithms,iterations);
counter = 0;
    
for i=1:iterations

    % generate experiment        
    [v1,v2,t,R] = create2D2DOmniExperiment(pt_number,cam_number,noise,outlier_fraction);
    [t_perturbed,R_perturbed] = perturb(t,R,0.01);
    T_perturbed = [R_perturbed,t_perturbed];
    T_init = [eye(3) zeros(3,1)];

    for a=1:num_algorithms
        tic
        Out = opengv_donotuse(algorithms{a},indices{a},v1,v2,T_init);
        execution_times(a,i) = toc/20.0;
    end

    counter = counter + 1;
    if counter == 1
        counter = 0;
        display(['Iteration ' num2str(i) ' of ' num2str(iterations) '(noise level ' num2str(noise) ')']);
    end        
end

%% Plot results

hist(log10(execution_times)')

legend(names,'Location','NorthWest')
xlabel('execution time [s]')
ylabel('occurence')
grid on

%% print the mean and median execution time on the console

display( 'mean execution times:' )
display(['sixpt:       ' num2str(mean(execution_times(1,:)'))] );
display(['ge:          ' num2str(mean(execution_times(2,:)'))] );
display(['seventeenpt: ' num2str(mean(execution_times(3,:)'))] );

%% Plot the results
% 
% [y1,x1] = hist(execution_times(1,:));
% [y2,x2] = hist(execution_times(2,:));
% [y3,x3] = hist(execution_times(3,:));
% 
% y1 = y1 / (x1(1,2) - x1(1,1));
% y2 = y2 / (x2(1,2) - x2(1,1));
% y3 = y3 / (x3(1,2) - x3(1,1));
% 
% figure(2)
% hold on
% plot(x1,y1,'b');
% plot(x2,y2,'g');
% plot(x3,y3,'r');
% legend(names,'Location','NorthWest')
% xlabel('execution time [s]')
% ylabel('probability')
% grid on
