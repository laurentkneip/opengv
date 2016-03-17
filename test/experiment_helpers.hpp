/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

#ifndef OPENGV_EXPERIMENT_HELPERS_HPP_
#define OPENGV_EXPERIMENT_HELPERS_HPP_

#include <opengv/types.hpp>
#include <memory>

namespace opengv
{

void generateCentralCameraSystem(
    translations_t & camOffsets,
    rotations_t & camRotations );

void generateRandomCameraSystem(
    int numberCameras,
    translations_t & camOffsets,
    rotations_t & camRotations );

void extractRelativePose(
    const translation_t & position1,
    const translation_t & position2,
    const rotation_t & rotation1,
    const rotation_t & rotation2,
    translation_t & relativePosition,
    rotation_t & relativeRotation,
    bool normalize = true );

void printExperimentCharacteristics(
    const translation_t & position,
    const rotation_t & rotation,
    double noise,
    double outlierFraction );

void printBearingVectorArraysMatlab(
    const bearingVectors_t & bearingVectors1,
    const bearingVectors_t & bearingVectors2 );

void printEssentialMatrix(
    const translation_t & position,
    const rotation_t & rotation);

void getPerturbedPose(
    const translation_t & position,
    const rotation_t & rotation,
    translation_t & perturbedPosition,
    rotation_t & perturbedRotation,
    double amplitude );

std::vector<int> getNindices( int n );




void generateRandom2D3DCorrespondences(
    const translation_t & position,
    const rotation_t & rotation,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    bearingVectors_t & bearingVectors,
    points_t & points,
    std::vector<int> & camCorrespondences,
    Eigen::MatrixXd & gt );

void generateMulti2D3DCorrespondences(
    const translation_t & position,
    const rotation_t & rotation,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t pointsPerCam,
    double noise,
    double outlierFraction,
    std::vector<std::shared_ptr<points_t> > & multiPoints,
    std::vector<std::shared_ptr<bearingVectors_t> > & multiBearingVectors,
    std::vector<std::shared_ptr<Eigen::MatrixXd> > & gt );

void generateRandom2D2DCorrespondences(
    const translation_t & position1,
    const rotation_t & rotation1,
    const translation_t & position2,
    const rotation_t & rotation2,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    bearingVectors_t & bearingVectors1,
    bearingVectors_t & bearingVectors2,
    std::vector<int> & camCorrespondences1,
    std::vector<int> & camCorrespondences2,
    Eigen::MatrixXd & gt );

void generateRandom3D3DCorrespondences(
    const translation_t & position1,
    const rotation_t & rotation1,
    const translation_t & position2,
    const rotation_t & rotation2,
    size_t numberPoints,
    double noise,
    double outlierFraction,
    bearingVectors_t & points1,
    bearingVectors_t & points2,
    Eigen::MatrixXd & gt );

void generateMulti2D2DCorrespondences(
    const translation_t & position1,
    const rotation_t & rotation1,
    const translation_t & position2,
    const rotation_t & rotation2,
    const translations_t & camOffsets,
    const rotations_t & camRotations,
    size_t pointsPerCam,
    double noise,
    double outlierFraction,
    std::vector<std::shared_ptr<bearingVectors_t> > & multiBearingVectors1,
    std::vector<std::shared_ptr<bearingVectors_t> > & multiBearingVectors2,
    std::vector<std::shared_ptr<Eigen::MatrixXd> > & gt );

}

#endif /* OPENGV_EXPERIMENT_HELPERS_HPP_ */
