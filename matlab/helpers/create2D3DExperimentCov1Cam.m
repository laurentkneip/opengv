function [points, v, cov, t, R ] = create2D3DExperimentCov1Cam(pt_number, noise, outlier_fraction)

%% Generate random point-cloud
minDepth = 4.0;
maxDepth = 8.0;

body_points  = [xrand(1,pt_number,[-2 2]); xrand(1,pt_number,[-2 2]); xrand(1,pt_number,[minDepth maxDepth])];
t   = mean(body_points,2);
R   = rodrigues(randn(3,1));
points = R\(body_points-repmat(t,1,pt_number));
            
%% Now create the correspondences by looping through the cameras

focal_length = 800.0;
K = [focal_length 0 0; 0 focal_length 0; 0 0 1];
cov = zeros(9,pt_number);

% project poitns
xx = [body_points(1,:)./body_points(3,:); body_points(2,:)./body_points(3,:)]*focal_length;
std = rand(1,pt_number)*noise;
noiseAdd = randn(2,pt_number).*[std;std];
xxn = xx+noiseAdd;
homx = [xxn/focal_length; ones(1,size(xxn,2))];
v = normc(homx);
for i=1:pt_number
    % covariance projection
    cov_proj = K\diag([std(i)^2 std(i)^2 0])/K';
    J = (eye(3)-(v(1:3,i)*v(1:3,i)')/(v(1:3,i)'*v(1:3,i)))/norm(homx);
    Evv = J*cov_proj*J';
    cov(:,i) = reshape(Evv,9,1);
end
%% copy over the position and orientation

R = R';
t = -R*t;