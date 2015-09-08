import pyopengv
import numpy as np


def generateRandomPoint( maximumDepth, minimumDepth ):
    cleanPoint = np.random.uniform(-1.0, 1.0)
    direction = cleanPoint / np.linalg.norm(cleanPoint)
    return (maximumDepth - minimumDepth) * cleanPoint + minimumDepth * direction


def generateRandomTranslation( maximumParallax ):
    return np.random.uniform(-maximumParallax, maximumParallax, 3)


def generateRandomRotation( maxAngle ):
    rpy = np.random.uniform(-maxAngle, maxAngle, 3)

    R1 = np.array([[1.0,  0.0,  0.0],
                   [0.0,  np.cos(rpy[0]), -np.sin(rpy[0])],
                   [0.0,  np.sin(rpy[0]),  np.cos(rpy[0])]])

    R2 = np.array([[ np.cos(rpy[1]),  0.0,  np.sin(rpy[1])],
                   [0.0,  1.0,  0.0],
                   [-np.sin(rpy[1]),  0.0,  np.cos(rpy[1])]])

    R3 = np.array([[np.cos(rpy[2]), -np.sin(rpy[2]),  0.0],
                   [np.sin(rpy[2]),  np.cos(rpy[2]),  0.0],
                   [0.0,  0.0,  1.0]])

    return R3.dot(R2.dot(R1))


def addNoise( noiseLevel, cleanPoint ):
    noisyPoint = cleanPoint + np.random.uniform(-noiseLevel, noiseLevel, 3)
    return noisyPoint / np.linalg.norm(noisyPoint)


def generateRandom2D2DCorrespondences(position1, rotation1,
                                      position2, rotation2,
                                      num_points,
                                      noise, outlier_fraction):
    min_depth = 4
    max_depth = 8

    # initialize point-cloud
    gt = np.empty((num_points, 3))
    for i in range(num_points):
        gt[i] = generateRandomPoint(max_depth, min_depth)

    bearing_vectors1 = np.empty((num_points, 3))
    bearing_vectors2 = np.empty((num_points, 3))
    for i in range(num_points):
        # get the point in viewpoint 1
        body_point = rotation1.T.dot(gt[i] - position1)

        # get the point in viewpoint 2
        body_point = rotation2.T.dot(gt[i] - position2)

        bearing_vectors1[i] = body_point / np.linalg.norm(body_point)
        bearing_vectors2[i] = body_point / np.linalg.norm(body_point)

        # add noise
        if noise > 0.0:
            bearing_vectors1[i] = addNoise(noise, bearing_vectors1[i])
            bearing_vectors2[i] = addNoise(noise, bearing_vectors2[i])

    # add outliers
    num_outliers = int(outlier_fraction * num_points)
    for i in range(num_outliers):
        # create random point
        p = generateRandomPoint(max_depth, min_depth)

        # project this point into viewpoint 2
        body_point = rotation2.T.dot(p - position2)

        # normalize the bearing vector
        bearing_vectors2[i] = body_point / np.linalg.norm(body_point)

    return bearing_vectors1, bearing_vectors2


def extractRelativePose(position1, position2, rotation1, rotation2):
    relativeRotation = rotation1.T.dot(rotation2)
    relativePosition = rotation1.T.dot(position2 - position1)
    return relativePosition, relativeRotation


def essentialMatrix(position, rotation):
  # E transforms vectors from vp 2 to 1: x_1^T * E * x_2 = 0
  # and E = (t)_skew*R
  t_skew = np.zeros((3, 3))
  t_skew[0,1] = -position[2]
  t_skew[0,2] = position[1]
  t_skew[1,0] = position[2]
  t_skew[1,2] = -position[0]
  t_skew[2,0] = -position[1]
  t_skew[2,1] = position[0]

  E = t_skew.dot(rotation)
  return E / np.linalg.norm(E)


def getPerturbedPose(position, rotation, amplitude):
    dp = generateRandomTranslation(amplitude)
    dR = generateRandomRotation(amplitude)
    return position + dp, rotation.dot(dR)


def test_relative_pose_ransac():

    # set experiment parameters
    noise = 0.0
    outlier_fraction = 0.0
    num_points = 10

    # generate a random pose for viewpoint 1
    position1 = np.zeros(3)
    rotation1 = np.eye(3)

    # generate a random pose for viewpoint 2
    position2 = generateRandomTranslation(2.0)
    rotation2 = generateRandomRotation(0.5)

    # derive correspondences based on random point-cloud
    bearing_vectors1, bearing_vectors2 = generateRandom2D2DCorrespondences(
        position1, rotation1, position2, rotation2,
        num_points, noise, outlier_fraction)

    # Extract the relative pose
    position, rotation = extractRelativePose(
        position1, position2, rotation1, rotation2)

    # print experiment characteristics
    print "the random position is:"
    print position
    print
    print "the random rotation is:"
    print rotation
    print
    print "the noise in the data is:"
    print noise
    print "the outlier fraction is:"
    print outlier_fraction

    # compute and print the essential-matrix
    print "the ground truth essential matrix is:"
    print essentialMatrix(position, rotation)
    print


    # running experiments
    print "running twopt"
    twopt_translation = pyopengv.relative_pose_twopt(bearing_vectors1, bearing_vectors2, rotation)

    print "running fivept_nister"
    fivept_nister_essentials = pyopengv.relative_pose_fivept_nister(bearing_vectors1, bearing_vectors2)

    print "running fivept_kneip"
    fivept_kneip_rotations = pyopengv.relative_pose_fivept_kneip(bearing_vectors1, bearing_vectors2)

    print "running sevenpt"
    sevenpt_essentials = pyopengv.relative_pose_sevenpt(bearing_vectors1, bearing_vectors2)

    print "running eightpt"
    eightpt_essential = pyopengv.relative_pose_eightpt(bearing_vectors1, bearing_vectors2)

    print "setting perturbed rotation and running eigensolver"
    t_perturbed, R_perturbed = getPerturbedPose( position, rotation, 0.01)
    eigensolver_rotation = pyopengv.relative_pose_eigensolver(bearing_vectors1, bearing_vectors2, R_perturbed)

    print "setting perturbed pose and performing nonlinear optimization"
    t_perturbed, R_perturbed = getPerturbedPose( position, rotation, 0.1)
    nonlinear_transformation = pyopengv.relative_pose_optimize_nonlinear(bearing_vectors1, bearing_vectors2, t_perturbed, R_perturbed)

    # print results
    print "results from two-points algorithm:"
    print twopt_translation
    print
    print "results from nisters' five-point algorithm:"
    for E in fivept_nister_essentials:
        print E
    print "results from kneip's five-point algorithm:"
    for R in fivept_kneip_rotations:
        print R
    print "results from seven-point algorithm:"
    for E in sevenpt_essentials:
        print E
    print "results from eight-point algorithm:"
    print eightpt_essential
    print "results from eigensystem based rotation solver:"
    print eigensolver_rotation
    print "results from nonlinear algorithm:"
    print nonlinear_transformation


test_relative_pose_ransac()
