%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

rng shuffle

%% Configure the benchmark

% noncentral case
cam_number = 4;
% Getting 10 points, and testing all algorithms with the respective number of points
pt_per_cam = 20;
% outlier test, so constant noise
noise = 0.5;
% repeat 100 tests per outlier-ratio
iterations = 50;

% The algorithms we want to test
algorithms = [ 0; 1; 2 ];
% The name of the algorithms in the final plots
names = { '6pt'; 'ge (8pt)'; '17pt'};

% The main experiment parameters
min_outlier_fraction = 0.05;
max_outlier_fraction = 0.45;
outlier_fraction_step = 0.05;

%% Run the benchmark

%prepare the overall result arrays
number_outlier_fraction_levels = round((max_outlier_fraction - min_outlier_fraction) / outlier_fraction_step + 1);
num_algorithms = size(algorithms,1);

%Run the experiment
for n=1:number_outlier_fraction_levels

    outlier_fraction = (n - 1) * outlier_fraction_step + min_outlier_fraction;
    display(['Analyzing outlier fraction level: ' num2str(outlier_fraction)])
    
    clear number_iterations
    clear execution_times
    
    counter = 0;
    
    temp_file_name1 = ['number_iterations_' num2str(outlier_fraction) '.mat'];
    temp_file_name2 = ['execution_times_' num2str(outlier_fraction) '.mat'];
    
    if exist(temp_file_name1,'file') > 0
        display(['number_iterations_' num2str(outlier_fraction) '.mat exists already'])
        load(temp_file_name1)
        load(temp_file_name2)
        startingIteration = size(number_iterations,2) + 1;
        display(['starting at ' num2str(startingIteration)])
    else
        startingIteration = 1;
    end
    
    if startingIteration <= iterations
        for i=startingIteration:iterations
        
            % generate experiment        
            [v1,v2,cam_offsets,t,R] = createMulti2D2DOmniExperiment(pt_per_cam,cam_number,noise,outlier_fraction);

            for a=1:num_algorithms
                
                if strcmp(names{a,1},'6pt') && outlier_fraction > 0.25
                    Out = zeros(4,5);
                    time = 10000000.0;
                else
                    tic
                    Out = opengv_experimental1( v1{1,1}, v1{2,1}, v1{3,1}, v1{4,1}, v2{1,1}, v2{2,1}, v2{3,1}, v2{4,1}, cam_offsets, algorithms(a,1) );
                    time = toc;
                end

                number_iterations(a,i) = Out(1,5);
                execution_times(a,i) = time;
            end

            save(temp_file_name1,'number_iterations');
            save(temp_file_name2,'execution_times');

            counter = counter + 1;
            if counter == 1
                counter = 0;
                display(['Iteration ' num2str(i) ' of ' num2str(iterations) '(outlier_fraction level ' num2str(outlier_fraction) ')']);
            end
        end
    end
end