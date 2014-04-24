%% Reset everything

clear all;
clc;
close all;
addpath('helpers');

%% Configure the benchmark

% noncentral case
cam_number = 4;
% Getting 20 points
pt_number = 20;
% no outliers
outlier_fraction = 0.0;
% no noise
noise = 0.0;

[v1,v2,t,R] = create2D2DOmniExperiment(pt_number,cam_number,noise,outlier_fraction);


%% Plot the smallest Eigenvalue

evs = zeros(101,101);

for cay_z_index=1:5
    
    cay_z = (cay_z_index - 3) * 0.2;
    cay_z_index
    
    for i=1:101
        cay_x = 0.4*((i-51)/50);
        i
        for j=1:101
        
            cay_y = 0.4*((j-51) / 50);
            
            x = [cay_x;cay_y;cay_z];
            
            if i == 51 && j == 51 && cay_z_index == 3
                cay_y = 0.4*((52-51) / 50);
                x = [cay_x;cay_y;cay_z];
            end
            
            Out = opengv('ge2',[1:1:size(v1,2)],v1,v2,[cayley2rot(x) zeros(3,1)]);
            evs(i,j) = Out(4,5);
            
        end
    end
    
    negativeIndices = evs<0;
    evs(negativeIndices) = 0.0;
    
    figure(1)
    hold on
    surf([-0.4 0.4],[-0.4 0.4],repmat(cay_z, [2 2]),log(evs),'facecolor','texture')
    %colormap(gray);
    view(45,30);
    daspect([2.5 2.5 1]);
    axis([-0.4 0.4 -0.4 0.4 -0.4 0.4])
    grid on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar

end

%% insert the zero layer
insertZeroLayer = 0;

if insertZeroLayer == 1
    
    gt = rot2cayley(R);
    cay_z = gt(3,1);
    
    for i=1:101
        cay_x = 0.4*((i-51)/50);
        i
        for j=1:101
        
            cay_y = 0.4*((j-51) / 50);
            
            x = [cay_x;cay_y;cay_z];
            
            if i == 51 && j == 51 && cay_z_index == 3
                cay_y = 0.4*((52-51) / 50);
                x = [cay_x;cay_y;cay_z];
            end
            
            Out = opengv('ge2',[1:1:size(v1,2)],v1,v2,[cayley2rot(x) zeros(3,1)]);
            evs(i,j) = Out(4,5);
            
        end
    end
    
    negativeIndices = evs<0;
    evs(negativeIndices) = 0.0;
    
    figure(1)
    hold on
    surf([-0.4 0.4],[-0.4 0.4],repmat(cay_z, [2 2]),log(evs),'facecolor','texture')
    %colormap(gray);
    view(45,30);
    daspect([2.5 2.5 1]);
    axis([-0.4 0.4 -0.4 0.4 -0.4 0.4])
    grid on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar
    
end

%% Plot the second smallest Eigenvalue

evs = zeros(101,101);

for cay_z_index=1:5
    
    cay_z = (cay_z_index - 3) * 0.2;
    cay_z_index
    
    for i=1:101
        cay_x = 0.4*((i-51)/50);
        i
        for j=1:101
        
            cay_y = 0.4*((j-51) / 50);
            
            x = [cay_x;cay_y;cay_z];
            
            if i == 51 && j == 51 && cay_z_index == 3
                cay_y = 0.4*((52-51) / 50);
                x = [cay_x;cay_y;cay_z];
            end
            
            Out = opengv('ge2',[1:1:size(v1,2)],v1,v2,[cayley2rot(x) zeros(3,1)]);
            evs(i,j) = Out(3,5);
            
        end
    end
    
    negativeIndices = evs<0;
    evs(negativeIndices) = 0.0;
    
    figure(2)
    hold on
    surf([-0.4 0.4],[-0.4 0.4],repmat(cay_z, [2 2]),log(evs),'facecolor','texture')
    %colormap(gray);
    view(45,30);
    daspect([2.5 2.5 1]);
    axis([-0.4 0.4 -0.4 0.4 -0.4 0.4])
    grid on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar

end

%% insert the zero layer
insertZeroLayer = 0;

if insertZeroLayer == 1
    
    gt = rot2cayley(R);
    cay_z = gt(3,1);
    
    for i=1:101
        cay_x = 0.4*((i-51)/50);
        i
        for j=1:101
        
            cay_y = 0.4*((j-51) / 50);
            
            x = [cay_x;cay_y;cay_z];
            
            if i == 51 && j == 51 && cay_z_index == 3
                cay_y = 0.4*((52-51) / 50);
                x = [cay_x;cay_y;cay_z];
            end
            
            Out = opengv('ge2',[1:1:size(v1,2)],v1,v2,[cayley2rot(x) zeros(3,1)]);
            evs(i,j) = Out(3,5);
            
        end
    end
    
    negativeIndices = evs<0;
    evs(negativeIndices) = 0.0;
    
    figure(2)
    hold on
    surf([-0.4 0.4],[-0.4 0.4],repmat(cay_z, [2 2]),log(evs),'facecolor','texture')
    %colormap(gray);
    view(45,30);
    daspect([2.5 2.5 1]);
    axis([-0.4 0.4 -0.4 0.4 -0.4 0.4])
    grid on
    xlabel('x');
    ylabel('y');
    zlabel('z');
    colorbar
    
end

%% Print the ground truth value on the console

minimum = rot2cayley(R)