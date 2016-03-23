function [points, v, t, R ] = create2D3DExperiment( pt_number, cam_number, noise, outlier_fraction )

%% generate the camera system

avg_cam_distance = 0.5;
cam_offsets = zeros(3,cam_number);
%cam_rotations = zeros(3,cam_number*3);

if cam_number == 1
    cam_offsets = zeros(3,1);
    %cam_rotations = eye(3);
else
    for i=1:cam_number
        cam_offsets(:,i) = avg_cam_distance * generateRandomR() * [1.0; 0.0; 0.0];
        %cam_rotations(:,(i-1)*3+1:(i-1)*3+3) = generateRandomR();
    end
end

%% generate random view-points

max_parallax = 2.0;
max_rotation = 0.5;

position = max_parallax * 2.0 * (rand(3,1) - repmat(0.5,3,1));
rotation = generateBoundedR(max_rotation);

%% Generate random point-cloud

minDepth = 4.0;
maxDepth = 8.0;

normalizedPoints = 2.0*(rand(3,pt_number)-repmat(0.5,3,pt_number));
norms = sqrt(sum(normalizedPoints.*normalizedPoints));
directions = normalizedPoints./repmat(norms,3,1);
points = (maxDepth-minDepth) * normalizedPoints + minDepth * directions;

%% Now create the correspondences by looping through the cameras

focal_length = 800.0;

v = zeros(6,pt_number);
cam_correspondence = 1;
cam_correspondences = zeros(1,pt_number);

for i=1:pt_number
    
    cam_offset = cam_offsets(:,cam_correspondence);
    %cam_rotation = cam_rotations(:,(cam_correspondence-1)*3+1:(cam_correspondence-1)*3+3);
    
    body_point = rotation' * (points(:,i)-position);
    
    % we actually omit the can rotation here by unrotating the bearing
    % vectors already
    bearingVector = body_point - cam_offset;
    bearingVector_norm = norm(bearingVector);
    bearingVector = bearingVector/bearingVector_norm;
    
    % add noise to the bearing vectors here
    bearingVector_noisy = addNoise(bearingVector,focal_length,noise);
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector_norm = norm(bearingVector_noisy);
    
    v(:,i) = [bearingVector_noisy./bearingVector_norm; cam_offset];
    
    % change the camera correspondence
    cam_correspondences(1,i) = cam_correspondence;
    cam_correspondence = cam_correspondence + 1;
    if cam_correspondence > cam_number
        cam_correspondence = 1;
    end
end

%% Add outliers
number_outliers = floor(outlier_fraction*pt_number);

if number_outliers > 0
for i=1:number_outliers
    
    cam_correspondence = cam_correspondences(1,i);
    
    cam_offset = cam_offsets(:,cam_correspondence);
    %cam_rotation = cam_rotations(:,(cam_correspondence-1)*3+1:(cam_correspondence-1)*3+3);
    
    %generate random point
    normalizedPoint = 2.0*(rand(3,1)-repmat(0.5,3,1));
    norm1 = sqrt(sum(normalizedPoint.*normalizedPoint));
    direction = normalizedPoint./norm1;
    point = (maxDepth-minDepth) * normalizedPoint + minDepth * direction;
    
    body_point = rotation' * (point-position);
    
    % store the point (no need to add noise)
    bearingVector = body_point - cam_offset;
    
    % store the normalized bearing vectors along with the cameras they are
    % being seen (we create correspondences that always originate from the
    % same camera, you can change this if you want)
    bearingVector_norm = norm(bearingVector);
    
    v(:,i) = [bearingVector./bearingVector_norm; cam_offset];
end
end

%% copy over the position and orientation

t = position;
R = rotation;

%% cut the cam offsets in the single camera (e.g. central case)

if cam_number == 1
    v = v(1:3,:);
end