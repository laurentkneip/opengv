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


#include <opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/triangulation/methods.hpp>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
#include <eigen3/unsupported/Eigen/NumericalDiff>

#include <opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp>

bool
opengv::sac_problems::
    relative_pose::NoncentralRelativePoseSacProblem::computeModelCoefficients(
    const std::vector<int> &indices,
    model_t & outModel) const
{
  bool returnValue = true;

  if(!_asCentral)
    outModel = opengv::relative_pose::seventeenpt(_adapter,indices);
  else
  {
    typedef opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem::algorithm_t algorithm_t;
    algorithm_t algorithm =
        opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem::STEWENIUS;
    opengv::sac_problems::relative_pose::CentralRelativePoseSacProblem centralProblem(_adapter, algorithm);
    returnValue = centralProblem.computeModelCoefficients(indices,outModel);

    //The transformation has been computed from cam to cam now, so transform that into the body frame
    translation_t t_c1c2 = outModel.col(3);
    rotation_t R_c1c2 = outModel.block<3,3>(0,0);
    translation_t t_bc1 = _adapter.getCamOffset1(indices[0]);
    rotation_t R_bc1 = _adapter.getCamRotation1(indices[0]);
    translation_t t_bc2 = _adapter.getCamOffset2(indices[0]);
    rotation_t R_bc2 = _adapter.getCamRotation2(indices[0]);

    outModel.block<3,3>(0,0) = R_bc1 * R_c1c2 * R_bc2.transpose();
    outModel.col(3) = t_bc1 + R_bc1 * t_c1c2 - outModel.block<3,3>(0,0) * t_bc2;
  }

  return returnValue;
}

void
opengv::sac_problems::
    relative_pose::NoncentralRelativePoseSacProblem::getSelectedDistancesToModel(
    const model_t & model,
    const std::vector<int> & indices,
    std::vector<double> & scores) const
{
  translation_t translation = model.col(3);
  rotation_t rotation = model.block<3,3>(0,0);

  Eigen::Matrix<double,4,1> p_hom;
  p_hom[3] = 1.0;

  for( size_t i = 0; i < indices.size(); i++ )
  {
    translation_t cam1Offset = _adapter.getCamOffset1(indices[i]);
    rotation_t cam1Rotation = _adapter.getCamRotation1(indices[i]);
    translation_t cam2Offset = _adapter.getCamOffset2(indices[i]);
    rotation_t cam2Rotation = _adapter.getCamRotation2(indices[i]);

    translation_t directTranslation =
        cam1Rotation.transpose() *
        ((translation - cam1Offset) + rotation * cam2Offset);
    rotation_t directRotation =
        cam1Rotation.transpose() * rotation * cam2Rotation;

    _adapter.sett12(directTranslation);
    _adapter.setR12(directRotation);

    transformation_t inverseSolution;
    inverseSolution.block<3,3>(0,0) = directRotation.transpose();
    inverseSolution.col(3) =
        -inverseSolution.block<3,3>(0,0)*directTranslation;

    p_hom.block<3,1>(0,0) =
        opengv::triangulation::triangulate2(_adapter,indices[i]);
    bearingVector_t reprojection1 = p_hom.block<3,1>(0,0);
    bearingVector_t reprojection2 = inverseSolution * p_hom;
    reprojection1 = reprojection1 / reprojection1.norm();
    reprojection2 = reprojection2 / reprojection2.norm();
    bearingVector_t f1 = _adapter.getBearingVector1(indices[i]);
    bearingVector_t f2 = _adapter.getBearingVector2(indices[i]);

    //bearing-vector based outlier criterium (select threshold accordingly):
    //1-(f1'*f2) = 1-cos(alpha) \in [0:2]
    double reprojError1 = 1.0 - (f1.transpose() * reprojection1);
    double reprojError2 = 1.0 - (f2.transpose() * reprojection2);
    scores.push_back(reprojError1 + reprojError2);
  }
}

void
opengv::sac_problems::
    relative_pose::NoncentralRelativePoseSacProblem::optimizeModelCoefficients(
    const std::vector<int> & inliers,
    const model_t & model,
    model_t & optimized_model)
{
  translation_t translation = model.col(3);
  rotation_t rotation = model.block<3,3>(0,0);
  _adapter.sett12(translation);
  _adapter.setR12(rotation);
  optimized_model =
      opengv::relative_pose::optimize_nonlinear(_adapter,inliers);
}

int
opengv::sac_problems::
    relative_pose::NoncentralRelativePoseSacProblem::getSampleSize() const
{
  int sampleSize = 17;
  if(_asCentral)
    sampleSize = 5 + 3;
  return sampleSize;
}
