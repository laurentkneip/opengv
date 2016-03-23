%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

% central case -> only one camera
cam_number = 1;
% let's only get 6 points, and generate new ones in each iteration
pt_number = 100;
% noise test, so no outliers
outlier_fraction = 0.1;

noise = 0.0;
[points,v,t,R] = create2D3DExperiment(pt_number,cam_number,noise,outlier_fraction);
[X, inliers] = opengv('p3p_kneip_ransac',points,v);