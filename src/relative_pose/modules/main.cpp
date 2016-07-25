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


#include <math.h>
#include <vector>
#include <Eigen/NonLinearOptimization>
#include <Eigen/NumericalDiff>

#include <opengv/relative_pose/modules/main.hpp>
#include <opengv/relative_pose/modules/fivept_nister/modules.hpp>
#include <opengv/relative_pose/modules/fivept_stewenius/modules.hpp>
#include <opengv/relative_pose/modules/fivept_kneip/modules.hpp>
#include <opengv/relative_pose/modules/eigensolver/modules.hpp>
#include <opengv/relative_pose/modules/sixpt/modules.hpp>
#include <opengv/relative_pose/modules/ge/modules.hpp>
#include <opengv/relative_pose/modules/sixpt_ventura/approx_relpose_generalized_fast_computeA.h>
#include <opengv/relative_pose/modules/sixpt_ventura/Polynomial.hpp>
#include <opengv/OptimizationFunctor.hpp>
#include <opengv/math/arun.hpp>
#include <opengv/math/cayley.hpp>
#include <opengv/math/Sturm.hpp>

#include <stdio.h>
#include <iostream>

void
opengv::relative_pose::modules::fivept_stewenius_main(
    const Eigen::Matrix<double,9,4> & EE,
    complexEssentials_t & complexEssentials )
{
  Eigen::Matrix<double,10,20> A;
  fivept_stewenius::composeA(EE,A);

  Eigen::Matrix<double,10,10> A1 = A.block(0,0,10,10);
  Eigen::Matrix<double,10,10> A2 = A.block(0,10,10,10);

  Eigen::FullPivLU< Eigen::Matrix<double,10,10> > luA1(A1);
  Eigen::Matrix<double,10,10> A3 = luA1.inverse()*A2;

  Eigen::Matrix<double,10,10> M;
  M.block(0,0,3,10) = -A3.block(0,0,3,10);
  M.block(3,0,2,10) = -A3.block(4,0,2,10);
  M.row(5) = -A3.row(7);
  M.block(6,0,4,10) = Eigen::Matrix<double,4,10>::Zero();
  M(6,0) = 1.0;
  M(7,1) = 1.0;
  M(8,3) = 1.0;
  M(9,6) = 1.0;

  Eigen::EigenSolver< Eigen::Matrix<double,10,10> > Eig(M,true);

  //Eigen::Matrix<std::complex<double>,10,1> D = Eig.eigenvalues();
  Eigen::Matrix<std::complex<double>,10,10> V = Eig.eigenvectors();

  Eigen::Matrix<std::complex<double>,3,10> V1;
  V1 = V.block(6,0,3,10);

  Eigen::Matrix<std::complex<double>,3,10> V2;
  V2.row(0) = V.row(9);
  V2.row(1) = V.row(9);
  V2.row(2) = V.row(9);

  Eigen::Matrix<std::complex<double>,4,10> SOLS;

  for( int r = 0; r < 3; r++ )
  {
    for( int c = 0; c < 10; c++ )
      SOLS(r,c) = V1(r,c)/V2(r,c);
  }

  SOLS.row(3) = Eigen::Matrix<std::complex<double>,1,10>::
      Constant(std::complex<double>(1.0,0.0));

  Eigen::Matrix<std::complex<double>,9,10> Evec = EE * SOLS;

  Eigen::Matrix<std::complex<double>,1,10> norms;
  for( int c = 0; c < 10; c++ )
  {
    norms(0,c) = std::complex<double>(0.0,0.0);

    for( int r = 0; r < 9; r++ )
      norms(0,c) += pow( Evec(r,c), 2 );

    norms(0,c) = sqrt(norms(0,c));
  }

  Eigen::Matrix<std::complex<double>,9,10> EvecNorms;
  for( int i = 0; i < 9; i++ )
    EvecNorms.row(i) = norms;

  for( int r = 0; r < 9; r++ )
  {
    for( int c = 0; c < 10; c++ )
      Evec(r,c) = Evec(r,c) / EvecNorms(r,c);
  }

  for( int c = 0; c < 10; c++ )
  {
    complexEssential_t complexEssential;

    complexEssential.row(0) = Evec.block(0,c,3,1).transpose();
    complexEssential.row(1) = Evec.block(3,c,3,1).transpose();
    complexEssential.row(2) = Evec.block(6,c,3,1).transpose();

    complexEssentials.push_back(complexEssential);
  }
}

void
opengv::relative_pose::modules::fivept_nister_main(
    const Eigen::Matrix<double,9,4> & EE,
    essentials_t & essentials )
{
  Eigen::Matrix<double,10,20> A;
  fivept_nister::composeA(EE,A);

  Eigen::Matrix<double,10,10> A1 = A.block(0,0,10,10);
  Eigen::Matrix<double,10,10> A2 = A.block(0,10,10,10);

  Eigen::FullPivLU< Eigen::Matrix<double,10,10> > luA1(A1);
  Eigen::Matrix<double,10,10> A3 = luA1.inverse()*A2;

  Eigen::Matrix<double,1,4> b11_part1 = Eigen::Matrix<double,1,4>::Zero();
  b11_part1.block<1,3>(0,1) = A3.block<1,3>(4,0);
  Eigen::Matrix<double,1,4> b11_part2 = Eigen::Matrix<double,1,4>::Zero();
  b11_part2.block<1,3>(0,0) = A3.block<1,3>(5,0);
  Eigen::Matrix<double,1,4> b11 = b11_part1 - b11_part2;

  Eigen::Matrix<double,1,4> b21_part1 = Eigen::Matrix<double,1,4>::Zero();
  b21_part1.block<1,3>(0,1) = A3.block<1,3>(6,0);
  Eigen::Matrix<double,1,4> b21_part2 = Eigen::Matrix<double,1,4>::Zero();
  b21_part2.block<1,3>(0,0) = A3.block<1,3>(7,0);
  Eigen::Matrix<double,1,4> b21 = b21_part1 - b21_part2;

  Eigen::Matrix<double,1,4> b31_part1 = Eigen::Matrix<double,1,4>::Zero();
  b31_part1.block<1,3>(0,1) = A3.block<1,3>(8,0);
  Eigen::Matrix<double,1,4> b31_part2 = Eigen::Matrix<double,1,4>::Zero();
  b31_part2.block<1,3>(0,0) = A3.block<1,3>(9,0);
  Eigen::Matrix<double,1,4> b31 = b31_part1 - b31_part2;

  Eigen::Matrix<double,1,4> b12_part1 = Eigen::Matrix<double,1,4>::Zero();
  b12_part1.block<1,3>(0,1) = A3.block<1,3>(4,3);
  Eigen::Matrix<double,1,4> b12_part2 = Eigen::Matrix<double,1,4>::Zero();
  b12_part2.block<1,3>(0,0) = A3.block<1,3>(5,3);
  Eigen::Matrix<double,1,4> b12 = b12_part1 - b12_part2;

  Eigen::Matrix<double,1,4> b22_part1 = Eigen::Matrix<double,1,4>::Zero();
  b22_part1.block<1,3>(0,1) = A3.block<1,3>(6,3);
  Eigen::Matrix<double,1,4> b22_part2 = Eigen::Matrix<double,1,4>::Zero();
  b22_part2.block<1,3>(0,0) = A3.block<1,3>(7,3);
  Eigen::Matrix<double,1,4> b22 = b22_part1 - b22_part2;

  Eigen::Matrix<double,1,4> b32_part1 = Eigen::Matrix<double,1,4>::Zero();
  b32_part1.block<1,3>(0,1) = A3.block<1,3>(8,3);
  Eigen::Matrix<double,1,4> b32_part2 = Eigen::Matrix<double,1,4>::Zero();
  b32_part2.block<1,3>(0,0) = A3.block<1,3>(9,3);
  Eigen::Matrix<double,1,4> b32 = b32_part1 - b32_part2;

  Eigen::Matrix<double,1,5> b13_part1 = Eigen::Matrix<double,1,5>::Zero();
  b13_part1.block<1,4>(0,1) = A3.block<1,4>(4,6);
  Eigen::Matrix<double,1,5> b13_part2 = Eigen::Matrix<double,1,5>::Zero();
  b13_part2.block<1,4>(0,0) = A3.block<1,4>(5,6);
  Eigen::Matrix<double,1,5> b13 = b13_part1 - b13_part2;

  Eigen::Matrix<double,1,5> b23_part1 = Eigen::Matrix<double,1,5>::Zero();
  b23_part1.block<1,4>(0,1) = A3.block<1,4>(6,6);
  Eigen::Matrix<double,1,5> b23_part2 = Eigen::Matrix<double,1,5>::Zero();
  b23_part2.block<1,4>(0,0) = A3.block<1,4>(7,6);
  Eigen::Matrix<double,1,5> b23 = b23_part1 - b23_part2;

  Eigen::Matrix<double,1,5> b33_part1 = Eigen::Matrix<double,1,5>::Zero();
  b33_part1.block<1,4>(0,1) = A3.block<1,4>(8,6);
  Eigen::Matrix<double,1,5> b33_part2 = Eigen::Matrix<double,1,5>::Zero();
  b33_part2.block<1,4>(0,0) = A3.block<1,4>(9,6);
  Eigen::Matrix<double,1,5> b33 = b33_part1 - b33_part2;

  Eigen::Matrix<double,1,8> p1_part1;
  fivept_nister::computeSeventhOrderPolynomial(b23,b12,p1_part1);
  Eigen::Matrix<double,1,8> p1_part2;
  fivept_nister::computeSeventhOrderPolynomial(b13,b22,p1_part2);
  Eigen::Matrix<double,1,8> p1 = p1_part1 - p1_part2;
  Eigen::Matrix<double,1,8> p2_part1;
  fivept_nister::computeSeventhOrderPolynomial(b13,b21,p2_part1);
  Eigen::Matrix<double,1,8> p2_part2;
  fivept_nister::computeSeventhOrderPolynomial(b23,b11,p2_part2);
  Eigen::Matrix<double,1,8> p2 = p2_part1 - p2_part2;
  Eigen::Matrix<double,1,7> p3_part1;
  fivept_nister::computeSixthOrderPolynomial(b11,b22,p3_part1);
  Eigen::Matrix<double,1,7> p3_part2;
  fivept_nister::computeSixthOrderPolynomial(b12,b21,p3_part2);
  Eigen::Matrix<double,1,7> p3 = p3_part1 - p3_part2;

  Eigen::Matrix<double,1,11> p_order10_part1;
  fivept_nister::computeTenthOrderPolynomialFrom73(p1,b31,p_order10_part1);
  Eigen::Matrix<double,1,11> p_order10_part2;
  fivept_nister::computeTenthOrderPolynomialFrom73(p2,b32,p_order10_part2);
  Eigen::Matrix<double,1,11> p_order10_part3;
  fivept_nister::computeTenthOrderPolynomialFrom64(p3,b33,p_order10_part3);
  Eigen::Matrix<double,1,11> p_order10 =
      p_order10_part1 + p_order10_part2 + p_order10_part3;

  math::Sturm sturmSequence(p_order10);
  std::vector<double> roots = sturmSequence.findRoots();

  Eigen::MatrixXd Evec(9,roots.size());

  for( size_t i = 0; i < roots.size(); i++ )
  {
    double z = roots[i];
    double x = fivept_nister::polyVal(p1,z)/fivept_nister::polyVal(p3,z);
    double y = fivept_nister::polyVal(p2,z)/fivept_nister::polyVal(p3,z);

    //pollishing here
    fivept_nister::pollishCoefficients(A,x,y,z);

    Evec.col(i) = x*EE.col(0) + y*EE.col(1) + z*EE.col(2) + EE.col(3);
  }

  Eigen::MatrixXd norms(1,roots.size());
  for( size_t c = 0; c < roots.size(); c++ )
  {
    norms(0,c) = 0.0;

    for( int r = 0; r < 9; r++ )
      norms(0,c) += pow( Evec(r,c), 2 );

    norms(0,c) = sqrt(norms(0,c));
  }

  Eigen::MatrixXd EvecNorms(9,roots.size());
  for( size_t i = 0; i < 9; i++ )
    EvecNorms.row(i) = norms;

  for( size_t r = 0; r < 9; r++ )
  {
    for( size_t c = 0; c < roots.size(); c++ )
      Evec(r,c) = Evec(r,c) / EvecNorms(r,c);
  }

  for( size_t c = 0; c < roots.size(); c++ )
  {
    essential_t essential;

    essential.row(0) = Evec.block<3,1>(0,c).transpose();
    essential.row(1) = Evec.block<3,1>(3,c).transpose();
    essential.row(2) = Evec.block<3,1>(6,c).transpose();

    essentials.push_back(essential);
  }
}

void
opengv::relative_pose::modules::fivept_kneip_main(
    const Eigen::Matrix<double,3,5> & f1,
    const Eigen::Matrix<double,3,5> & f2,
    rotations_t & rotations )
{
  Eigen::Matrix<double,66,197> groebnerMatrix =
      Eigen::Matrix<double,66,197>::Zero();
  Eigen::Matrix3d temp1 = Eigen::Matrix3d::Zero();
  Eigen::Matrix3d temp2 = Eigen::Matrix3d::Zero();
  std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> > c_1;
  c_1.push_back(temp1);
  c_1.push_back(temp1);
  c_1.push_back(temp1);
  std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> > c_2;
  c_2.push_back(temp2);
  c_2.push_back(temp2);
  c_2.push_back(temp2);

  int currentRow = 0;

  for( int firstFeat = 0; firstFeat < 5; firstFeat++ )
  {
    for( int secondFeat = firstFeat + 1; secondFeat < 5; secondFeat++ )
    {
      temp1 = f1.col(firstFeat)*f1.col(secondFeat).transpose();
      temp2 = f2.col(firstFeat)*f2.col(secondFeat).transpose();

      for( int thirdFeat = secondFeat + 1; thirdFeat < 5; thirdFeat++ )
      {
        c_1[0] = temp1 * f1(0,thirdFeat);
        c_1[1] = temp1 * f1(1,thirdFeat);
        c_1[2] = temp1 * f1(2,thirdFeat);
        c_2[0] = temp2 * f2(0,thirdFeat);
        c_2[1] = temp2 * f2(1,thirdFeat);
        c_2[2] = temp2 * f2(2,thirdFeat);

        groebnerMatrix.row(currentRow++) =
            fivept_kneip::initEpncpRowR( c_1, c_2 );
      }
    }
  }

  fivept_kneip::initMatrix(groebnerMatrix);
  fivept_kneip::computeBasis(groebnerMatrix);

  Eigen::Matrix<double,20,20> M = Eigen::Matrix<double,20,20>::Zero();
  M.block<10,20>(0,0) = -groebnerMatrix.block<10,20>(51,177);
  M(10,1) = 1.0;
  M(11,2) = 1.0;
  M(12,3) = 1.0;
  M(13,4) = 1.0;
  M(14,5) = 1.0;
  M(15,6) = 1.0;
  M(16,7) = 1.0;
  M(17,8) = 1.0;
  M(18,9) = 1.0;
  M(19,18) = 1.0;

  Eigen::EigenSolver< Eigen::Matrix<double,20,20> > Eig(M,true);
  Eigen::Matrix<std::complex<double>,20,1> D = Eig.eigenvalues();
  Eigen::Matrix<std::complex<double>,20,20> V = Eig.eigenvectors();

  //Eigen::Matrix<std::complex<double>,9,1> tempVector;
  rotation_t finalRotation = Eigen::Matrix3d::Zero();
  for( unsigned int i = 0; i < 20; i++ )
  {
    std::complex<double> tempp;
    tempp = D[i];

    //check if we have a real solution
    if( fabs(tempp.imag()) < 0.1 )
    {
      tempp = V(18,i)/V(19,i);
      finalRotation(0,0) = tempp.real();// tempVector[0] = tempp;
      tempp = V(17,i)/V(19,i);
      finalRotation(0,1) = tempp.real();// tempVector[1] = tempp;
      tempp = V(16,i)/V(19,i);
      finalRotation(0,2) = tempp.real();// tempVector[2] = tempp;
      tempp = V(15,i)/V(19,i);
      finalRotation(1,0) = tempp.real();// tempVector[3] = tempp;
      tempp = V(14,i)/V(19,i);
      finalRotation(1,1) = tempp.real();// tempVector[4] = tempp;
      tempp = V(13,i)/V(19,i);
      finalRotation(1,2) = tempp.real();// tempVector[5] = tempp;
      tempp = V(12,i)/V(19,i);
      finalRotation(2,0) = tempp.real();// tempVector[6] = tempp;
      tempp = V(11,i)/V(19,i);
      finalRotation(2,1) = tempp.real();// tempVector[7] = tempp;
      tempp = V(10,i)/V(19,i);
      finalRotation(2,2) = tempp.real();// tempVector[8] = tempp;

      double tempNorm = finalRotation.row(1).norm();
      finalRotation.row(1) = finalRotation.row(1) / tempNorm;
      tempNorm = finalRotation.row(2).norm();
      finalRotation.row(2) = finalRotation.row(2) / tempNorm;

      //check if the normalized rotation matrix has determinant close enough to 1
      if( fabs( finalRotation.determinant() - 1.0 ) < 0.1 )
      {
        Eigen::Matrix3d eval;
        double totalEval = 0.0;
        eval.col(0) = f1.col(0).cross(finalRotation*f2.col(0));
        eval.col(1) = f1.col(1).cross(finalRotation*f2.col(1));
        eval.col(2) = f1.col(2).cross(finalRotation*f2.col(2));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(0).cross(finalRotation*f2.col(0));
        eval.col(1) = f1.col(1).cross(finalRotation*f2.col(1));
        eval.col(2) = f1.col(3).cross(finalRotation*f2.col(3));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(0).cross(finalRotation*f2.col(0));
        eval.col(1) = f1.col(1).cross(finalRotation*f2.col(1));
        eval.col(2) = f1.col(4).cross(finalRotation*f2.col(4));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(0).cross(finalRotation*f2.col(0));
        eval.col(1) = f1.col(2).cross(finalRotation*f2.col(2));
        eval.col(2) = f1.col(3).cross(finalRotation*f2.col(3));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(0).cross(finalRotation*f2.col(0));
        eval.col(1) = f1.col(2).cross(finalRotation*f2.col(2));
        eval.col(2) = f1.col(4).cross(finalRotation*f2.col(4));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(0).cross(finalRotation*f2.col(0));
        eval.col(1) = f1.col(3).cross(finalRotation*f2.col(3));
        eval.col(2) = f1.col(4).cross(finalRotation*f2.col(4));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(1).cross(finalRotation*f2.col(1));
        eval.col(1) = f1.col(2).cross(finalRotation*f2.col(2));
        eval.col(2) = f1.col(3).cross(finalRotation*f2.col(3));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(1).cross(finalRotation*f2.col(1));
        eval.col(1) = f1.col(2).cross(finalRotation*f2.col(2));
        eval.col(2) = f1.col(4).cross(finalRotation*f2.col(4));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(1).cross(finalRotation*f2.col(1));
        eval.col(1) = f1.col(3).cross(finalRotation*f2.col(3));
        eval.col(2) = f1.col(4).cross(finalRotation*f2.col(4));
        totalEval += fabs(eval.determinant());
        eval.col(0) = f1.col(2).cross(finalRotation*f2.col(2));
        eval.col(1) = f1.col(3).cross(finalRotation*f2.col(3));
        eval.col(2) = f1.col(4).cross(finalRotation*f2.col(4));
        totalEval += fabs(eval.determinant());
        totalEval += fabs(
            finalRotation(0,0)*finalRotation(0,0)+
            finalRotation(0,1)*finalRotation(0,1)+
            finalRotation(0,2)*finalRotation(0,2)-1);
        totalEval += fabs(
            finalRotation(0,0)*finalRotation(1,0)+
            finalRotation(0,1)*finalRotation(1,1)+
            finalRotation(0,2)*finalRotation(1,2));
        totalEval += fabs(
            finalRotation(0,0)*finalRotation(2,0)+
            finalRotation(0,1)*finalRotation(2,1)+
            finalRotation(0,2)*finalRotation(2,2));
        totalEval += fabs(
            finalRotation(1,0)*finalRotation(1,0)+
            finalRotation(1,1)*finalRotation(1,1)+
            finalRotation(1,2)*finalRotation(1,2)-1);
        totalEval += fabs(
            finalRotation(1,0)*finalRotation(2,0)+
            finalRotation(1,1)*finalRotation(2,1)+
            finalRotation(1,2)*finalRotation(2,2));
        totalEval += fabs(
            finalRotation(2,0)*finalRotation(2,0)+
            finalRotation(2,1)*finalRotation(2,1)+
            finalRotation(2,2)*finalRotation(2,2)-1);
        totalEval += fabs(
            finalRotation(1,0)*finalRotation(2,1)-
            finalRotation(2,0)*finalRotation(1,1)-finalRotation(0,2));
        totalEval += fabs(
            finalRotation(2,0)*finalRotation(0,1)-
            finalRotation(0,0)*finalRotation(2,1)-finalRotation(1,2));
        totalEval += fabs(
            finalRotation(0,0)*finalRotation(1,1)-
            finalRotation(1,0)*finalRotation(0,1)-finalRotation(2,2));

        //check if the initial constraints are fullfilled to a sufficient extend
        if( totalEval < 0.001 )
        {
          Eigen::Matrix<double,3,5> normalVectors;
          normalVectors.col(0) = f1.col(0).cross(finalRotation * f2.col(0));
          normalVectors.col(1) = f1.col(1).cross(finalRotation * f2.col(1));
          normalVectors.col(2) = f1.col(2).cross(finalRotation * f2.col(2));
          normalVectors.col(3) = f1.col(3).cross(finalRotation * f2.col(3));
          normalVectors.col(4) = f1.col(4).cross(finalRotation * f2.col(4));

          Eigen::Vector3d trans01 =
              normalVectors.col(0).cross(normalVectors.col(1));
          Eigen::Vector3d trans02 =
              normalVectors.col(0).cross(normalVectors.col(2));
          Eigen::Vector3d trans03 =
              normalVectors.col(0).cross(normalVectors.col(3));
          Eigen::Vector3d trans04 =
              normalVectors.col(0).cross(normalVectors.col(4));
          Eigen::Vector3d trans12 =
              normalVectors.col(1).cross(normalVectors.col(2));
          Eigen::Vector3d trans13 =
              normalVectors.col(1).cross(normalVectors.col(3));
          Eigen::Vector3d trans14 =
              normalVectors.col(1).cross(normalVectors.col(4));
          Eigen::Vector3d trans23 =
              normalVectors.col(2).cross(normalVectors.col(3));
          Eigen::Vector3d trans24 =
              normalVectors.col(2).cross(normalVectors.col(4));
          Eigen::Vector3d trans34 =
              normalVectors.col(3).cross(normalVectors.col(4));

          Eigen::Vector3d tempVector1;
          Eigen::Vector3d tempVector2;

          bool positive = true;

          for( int i = 0; i < 5; i++ )
          {
            tempVector1 = trans01.cross( f1.col(i) );
            tempVector2 = trans01.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans02.cross( f1.col(i) );
            tempVector2 = trans02.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans03.cross( f1.col(i) );
            tempVector2 = trans03.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans04.cross( f1.col(i) );
            tempVector2 = trans04.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans12.cross( f1.col(i) );
            tempVector2 = trans12.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans13.cross( f1.col(i) );
            tempVector2 = trans13.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans14.cross( f1.col(i) );
            tempVector2 = trans14.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans23.cross( f1.col(i) );
            tempVector2 = trans23.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans24.cross( f1.col(i) );
            tempVector2 = trans24.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
            tempVector1 = trans34.cross( f1.col(i) );
            tempVector2 = trans34.cross( finalRotation * f2.col(i) );
            if( tempVector1.dot(tempVector2) < 0 )
            {
              positive = false;
              break;
            }
          }

          //finally, check if the cheiriality constraint is fullfilled to
          //sufficient extend
          if( positive )
            rotations.push_back(finalRotation);
        }
      }
    }
  }
}

namespace opengv
{
namespace relative_pose
{
namespace modules
{

struct Eigensolver_step : OptimizationFunctor<double>
{
  const Eigen::Matrix3d & _xxF;
  const Eigen::Matrix3d & _yyF;
  const Eigen::Matrix3d & _zzF;
  const Eigen::Matrix3d & _xyF;
  const Eigen::Matrix3d & _yzF;
  const Eigen::Matrix3d & _zxF;

  Eigensolver_step(
    const Eigen::Matrix3d & xxF,
    const Eigen::Matrix3d & yyF,
    const Eigen::Matrix3d & zzF,
    const Eigen::Matrix3d & xyF,
    const Eigen::Matrix3d & yzF,
    const Eigen::Matrix3d & zxF ) :
    OptimizationFunctor<double>(3,3),
    _xxF(xxF),_yyF(yyF),_zzF(zzF),_xyF(xyF),_yzF(yzF),_zxF(zxF) {}

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {
    cayley_t cayley = x;
    Eigen::Matrix<double,1,3> jacobian;
    eigensolver::getSmallestEVwithJacobian(
        _xxF,_yyF,_zzF,_xyF,_yzF,_zxF,cayley,jacobian);

    fvec[0] = jacobian(0,0);
    fvec[1] = jacobian(0,1);
    fvec[2] = jacobian(0,2);
    return 0;
  }
};

}
}
}

void
opengv::relative_pose::modules::eigensolver_main(
  const Eigen::Matrix3d & xxF,
  const Eigen::Matrix3d & yyF,
  const Eigen::Matrix3d & zzF,
  const Eigen::Matrix3d & xyF,
  const Eigen::Matrix3d & yzF,
  const Eigen::Matrix3d & zxF,
  eigensolverOutput_t & output )
{
  const int n=3;
  VectorXd x(n);

  x = math::rot2cayley(output.rotation);
  Eigensolver_step functor(xxF,yyF,zzF,xyF,yzF,zxF );
  NumericalDiff<Eigensolver_step> numDiff(functor);
  LevenbergMarquardt< NumericalDiff<Eigensolver_step> > lm(numDiff);

  lm.resetParameters();
  lm.parameters.ftol = 0.00005;
  lm.parameters.xtol = 1.E1*NumTraits<double>::epsilon();
  lm.parameters.maxfev = 100;
  lm.minimize(x);

  cayley_t cayley = x;
  rotation_t R = math::cayley2rot(cayley);

  Eigen::Matrix3d M = eigensolver::composeM(xxF,yyF,zzF,xyF,yzF,zxF,cayley);
  Eigen::EigenSolver< Eigen::Matrix3d > Eig(M,true);
  Eigen::Matrix<std::complex<double>,3,1> D_complex = Eig.eigenvalues();
  Eigen::Matrix<std::complex<double>,3,3> V_complex = Eig.eigenvectors();
  eigenvalues_t D;
  eigenvectors_t V;
  for(size_t i = 0; i < 3; i++)
  {
    D[i] = D_complex[i].real();
    for(size_t j = 0; j < 3; j++)
      V(i,j) = V_complex(i,j).real();
  }

  double translationMagnitude = sqrt(pow(D[1],2) + pow(D[2],2));
  translation_t t = translationMagnitude * V.col(0);

  output.translation = t;
  output.rotation = R;
  output.eigenvalues = D;
  output.eigenvectors = V;

}

void
opengv::relative_pose::modules::sixpt_main(
  Eigen::Matrix<double,6,6> & L1,
  Eigen::Matrix<double,6,6> & L2,
  rotations_t & solutions)
{
  
  //create vectors of Pluecker coordinates
  std::vector< Eigen::Matrix<double,6,1>, Eigen::aligned_allocator<Eigen::Matrix<double,6,1> > > L1vec;
  std::vector< Eigen::Matrix<double,6,1>, Eigen::aligned_allocator<Eigen::Matrix<double,6,1> > > L2vec;
  for(int i = 0; i < 6; i++ )
  {
    L1vec.push_back( L1.col(i) );
    L2vec.push_back( L2.col(i) );
  }
  
  //setup the action matrix
  Eigen::Matrix<double,64,64> Action = Eigen::Matrix<double,64,64>::Zero();
  sixpt::setupAction( L1vec, L2vec, Action );
  
  //finally eigen-decompose the action-matrix and obtain the solutions
  Eigen::EigenSolver< Eigen::Matrix<double,64,64> > Eig(Action,true);
  Eigen::Matrix<std::complex<double>,64,64> EV = Eig.eigenvectors();
  
  solutions.reserve(64);  
  for( int c = 0; c < 64; c++ )
  {
    cayley_t solution;
    for( int r = 0; r < 3; r++ )
    {
      std::complex<double> temp = EV(60+r,c)/EV(63,c);
      solution[r] = temp.real();
    }
    
    solutions.push_back(math::cayley2rot(solution).transpose());
  }
}

void
opengv::relative_pose::modules::sixpt_ventura_main(
const Eigen::Matrix<double, 6, 6> &w1,
const Eigen::Matrix<double, 6, 6> &w2,
const Eigen::Matrix<double, 6, 6> &w3,
const Eigen::Matrix<double, 6, 6> &w4,
const Eigen::Matrix<double, 6, 6> &w5,
const Eigen::Matrix<double, 6, 6> &w6,
std::vector<Eigen::Vector3d> &rsolns)
{
	const Eigen::Matrix<double, 15, 35> A = 
		opengv::relative_pose::modules::sixpt_ventura::computeA(w1, w2, w3, w4, w5, w6);

	const Eigen::Matrix<double, 15, 20> G = A.block<15, 15>(0, 0).lu().solve(A.block<15, 20>(0, 15));

	const double B11coeffs[] = { -G(6, 0), G(5, 0) - G(6, 1), G(5, 1) - G(6, 2), G(5, 2) };
	const double B12coeffs[] = { -G(6, 3), G(5, 3) - G(6, 4), G(5, 4) - G(6, 5), G(5, 5) };
	const double B13coeffs[] = { -G(6, 6), G(5, 6) - G(6, 7), G(5, 7) - G(6, 8), G(5, 8) };
	const double B14coeffs[] = { -G(6, 9), G(5, 9) - G(6, 10), G(5, 10) - G(6, 11), G(5, 11) - G(6, 12), G(5, 12) };
	const double B15coeffs[] = { -G(6, 13), G(5, 13) - G(6, 14), G(5, 14) - G(6, 15), G(5, 15) - G(6, 16), G(5, 16) };
	const double B16coeffs[] = { -G(6, 17), G(5, 17) - G(6, 18), G(5, 18) - G(6, 19), G(5, 19) };
	const double B21coeffs[] = { -G(8, 0), G(7, 0) - G(8, 1), G(7, 1) - G(8, 2), G(7, 2) };
	const double B22coeffs[] = { -G(8, 3), G(7, 3) - G(8, 4), G(7, 4) - G(8, 5), G(7, 5) };
	const double B23coeffs[] = { -G(8, 6), G(7, 6) - G(8, 7), G(7, 7) - G(8, 8), G(7, 8) };
	const double B24coeffs[] = { -G(8, 9), G(7, 9) - G(8, 10), G(7, 10) - G(8, 11), G(7, 11) - G(8, 12), G(7, 12) };
	const double B25coeffs[] = { -G(8, 13), G(7, 13) - G(8, 14), G(7, 14) - G(8, 15), G(7, 15) - G(8, 16), G(7, 16) };
	const double B26coeffs[] = { -G(8, 17), G(7, 17) - G(8, 18), G(7, 18) - G(8, 19), G(7, 19) };
	const double B31coeffs[] = { -G(10, 0), G(9, 0) - G(10, 1), G(9, 1) - G(10, 2), G(9, 2) };
	const double B32coeffs[] = { -G(10, 3), G(9, 3) - G(10, 4), G(9, 4) - G(10, 5), G(9, 5) };
	const double B33coeffs[] = { -G(10, 6), G(9, 6) - G(10, 7), G(9, 7) - G(10, 8), G(9, 8) };
	const double B34coeffs[] = { -G(10, 9), G(9, 9) - G(10, 10), G(9, 10) - G(10, 11), G(9, 11) - G(10, 12), G(9, 12) };
	const double B35coeffs[] = { -G(10, 13), G(9, 13) - G(10, 14), G(9, 14) - G(10, 15), G(9, 15) - G(10, 16), G(9, 16) };
	const double B36coeffs[] = { -G(10, 17), G(9, 17) - G(10, 18), G(9, 18) - G(10, 19), G(9, 19) };
	const double B41coeffs[] = { -G(12, 0), G(11, 0) - G(12, 1), G(11, 1) - G(12, 2), G(11, 2) };
	const double B42coeffs[] = { -G(12, 3), G(11, 3) - G(12, 4), G(11, 4) - G(12, 5), G(11, 5) };
	const double B43coeffs[] = { -G(12, 6), G(11, 6) - G(12, 7), G(11, 7) - G(12, 8), G(11, 8) };
	const double B44coeffs[] = { -G(12, 9), G(11, 9) - G(12, 10), G(11, 10) - G(12, 11), G(11, 11) - G(12, 12), G(11, 12) };
	const double B45coeffs[] = { -G(12, 13), G(11, 13) - G(12, 14), G(11, 14) - G(12, 15), G(11, 15) - G(12, 16), G(11, 16) };
	const double B46coeffs[] = { -G(12, 17), G(11, 17) - G(12, 18), G(11, 18) - G(12, 19), G(11, 19) };
	const double B51coeffs[] = { G(13, 0), G(13, 1), G(13, 2) };
	const double B52coeffs[] = { G(13, 3), G(13, 4), G(13, 5) };
	const double B53coeffs[] = { G(13, 6), G(13, 7), G(13, 8) };
	const double B54coeffs[] = { G(13, 9), G(13, 10), G(13, 11), G(13, 12) };
	const double B55coeffs[] = { G(13, 13), G(13, 14), G(13, 15), G(13, 16) };
	const double B56coeffs[] = { 1, 0, G(13, 17), G(13, 18), G(13, 19) };
	const double B61coeffs[] = { G(14, 0), G(14, 1), G(14, 2) };
	const double B62coeffs[] = { G(14, 3), G(14, 4), G(14, 5) };
	const double B63coeffs[] = { G(14, 6), G(14, 7), G(14, 8) };
	const double B64coeffs[] = { G(14, 9), G(14, 10), G(14, 11), G(14, 12) };
	const double B65coeffs[] = { G(14, 13), G(14, 14), G(14, 15), G(14, 16) };
	const double B66coeffs[] = { 1, G(14, 17), G(14, 18), G(14, 19) };
	const PolynomialVentura::Polynomial<3> B11(B11coeffs);
	const PolynomialVentura::Polynomial<3> B12(B12coeffs);
	const PolynomialVentura::Polynomial<3> B13(B13coeffs);
	const PolynomialVentura::Polynomial<4> B14(B14coeffs);
	const PolynomialVentura::Polynomial<4> B15(B15coeffs);
	const PolynomialVentura::Polynomial<3> B16(B16coeffs);
	const PolynomialVentura::Polynomial<3> B21(B21coeffs);
	const PolynomialVentura::Polynomial<3> B22(B22coeffs);
	const PolynomialVentura::Polynomial<3> B23(B23coeffs);
	const PolynomialVentura::Polynomial<4> B24(B24coeffs);
	const PolynomialVentura::Polynomial<4> B25(B25coeffs);
	const PolynomialVentura::Polynomial<3> B26(B26coeffs);
	const PolynomialVentura::Polynomial<3> B31(B31coeffs);
	const PolynomialVentura::Polynomial<3> B32(B32coeffs);
	const PolynomialVentura::Polynomial<3> B33(B33coeffs);
	const PolynomialVentura::Polynomial<4> B34(B34coeffs);
	const PolynomialVentura::Polynomial<4> B35(B35coeffs);
	const PolynomialVentura::Polynomial<3> B36(B36coeffs);
	const PolynomialVentura::Polynomial<3> B41(B41coeffs);
	const PolynomialVentura::Polynomial<3> B42(B42coeffs);
	const PolynomialVentura::Polynomial<3> B43(B43coeffs);
	const PolynomialVentura::Polynomial<4> B44(B44coeffs);
	const PolynomialVentura::Polynomial<4> B45(B45coeffs);
	const PolynomialVentura::Polynomial<3> B46(B46coeffs);
	const PolynomialVentura::Polynomial<2> B51(B51coeffs);
	const PolynomialVentura::Polynomial<2> B52(B52coeffs);
	const PolynomialVentura::Polynomial<2> B53(B53coeffs);
	const PolynomialVentura::Polynomial<3> B54(B54coeffs);
	const PolynomialVentura::Polynomial<3> B55(B55coeffs);
	const PolynomialVentura::Polynomial<4> B56(B56coeffs);
	const PolynomialVentura::Polynomial<2> B61(B61coeffs);
	const PolynomialVentura::Polynomial<2> B62(B62coeffs);
	const PolynomialVentura::Polynomial<2> B63(B63coeffs);
	const PolynomialVentura::Polynomial<3> B64(B64coeffs);
	const PolynomialVentura::Polynomial<3> B65(B65coeffs);
	const PolynomialVentura::Polynomial<3> B66(B66coeffs);
	const PolynomialVentura::Polynomial<6> t2 = B11*B22;
	const PolynomialVentura::Polynomial<6> t14 = B12*B21;
	const PolynomialVentura::Polynomial<6> t3 = t2 - t14;
	const PolynomialVentura::Polynomial<6> t4 = B11*B32;
	const PolynomialVentura::Polynomial<6> t16 = B12*B31;
	const PolynomialVentura::Polynomial<6> t5 = t4 - t16;
	const PolynomialVentura::Polynomial<6> t6 = B11*B42;
	const PolynomialVentura::Polynomial<6> t27 = B12*B41;
	const PolynomialVentura::Polynomial<6> t7 = t6 - t27;
	const PolynomialVentura::Polynomial<6> t8 = B21*B32;
	const PolynomialVentura::Polynomial<6> t17 = B22*B31;
	const PolynomialVentura::Polynomial<6> t9 = t8 - t17;
	const PolynomialVentura::Polynomial<6> t10 = B21*B42;
	const PolynomialVentura::Polynomial<6> t28 = B22*B41;
	const PolynomialVentura::Polynomial<6> t11 = t10 - t28;
	const PolynomialVentura::Polynomial<6> t12 = B31*B42;
	const PolynomialVentura::Polynomial<6> t39 = B32*B41;
	const PolynomialVentura::Polynomial<6> t13 = t12 - t39;
	const PolynomialVentura::Polynomial<9> t15 = B33*t3;
	const PolynomialVentura::Polynomial<9> t18 = B13*t9;
	const PolynomialVentura::Polynomial<9> t62 = B23*t5;
	const PolynomialVentura::Polynomial<9> t19 = t15 + t18 - t62;
	const PolynomialVentura::Polynomial<5> t20 = B11*B52;
	const PolynomialVentura::Polynomial<5> t32 = B12*B51;
	const PolynomialVentura::Polynomial<5> t21 = t20 - t32;
	const PolynomialVentura::Polynomial<5> t22 = B21*B52;
	const PolynomialVentura::Polynomial<5> t33 = B22*B51;
	const PolynomialVentura::Polynomial<5> t23 = t22 - t33;
	const PolynomialVentura::Polynomial<5> t24 = B31*B52;
	const PolynomialVentura::Polynomial<5> t43 = B32*B51;
	const PolynomialVentura::Polynomial<5> t25 = t24 - t43;
	const PolynomialVentura::Polynomial<9> t26 = B43*t3;
	const PolynomialVentura::Polynomial<9> t29 = B13*t11;
	const PolynomialVentura::Polynomial<9> t129 = B23*t7;
	const PolynomialVentura::Polynomial<9> t30 = t26 + t29 - t129;
	const PolynomialVentura::Polynomial<8> t31 = B53*t3;
	const PolynomialVentura::Polynomial<8> t34 = B13*t23;
	const PolynomialVentura::Polynomial<8> t63 = B23*t21;
	const PolynomialVentura::Polynomial<8> t35 = t31 + t34 - t63;
	const PolynomialVentura::Polynomial<5> t36 = B41*B52;
	const PolynomialVentura::Polynomial<5> t47 = B42*B51;
	const PolynomialVentura::Polynomial<5> t37 = t36 - t47;
	const PolynomialVentura::Polynomial<9> t38 = B43*t5;
	const PolynomialVentura::Polynomial<9> t40 = B13*t13;
	const PolynomialVentura::Polynomial<9> t99 = B33*t7;
	const PolynomialVentura::Polynomial<9> t41 = t38 + t40 - t99;
	const PolynomialVentura::Polynomial<8> t42 = B53*t5;
	const PolynomialVentura::Polynomial<8> t44 = B13*t25;
	const PolynomialVentura::Polynomial<8> t65 = B33*t21;
	const PolynomialVentura::Polynomial<8> t45 = t42 + t44 - t65;
	const PolynomialVentura::Polynomial<8> t46 = B53*t7;
	const PolynomialVentura::Polynomial<8> t48 = B13*t37;
	const PolynomialVentura::Polynomial<8> t101 = B43*t21;
	const PolynomialVentura::Polynomial<8> t49 = t46 + t48 - t101;
	const PolynomialVentura::Polynomial<9> t50 = B43*t9;
	const PolynomialVentura::Polynomial<9> t51 = B23*t13;
	const PolynomialVentura::Polynomial<9> t131 = B33*t11;
	const PolynomialVentura::Polynomial<9> t52 = t50 + t51 - t131;
	const PolynomialVentura::Polynomial<8> t53 = B53*t9;
	const PolynomialVentura::Polynomial<8> t54 = B23*t25;
	const PolynomialVentura::Polynomial<8> t66 = B33*t23;
	const PolynomialVentura::Polynomial<8> t55 = t53 + t54 - t66;
	const PolynomialVentura::Polynomial<8> t56 = B53*t11;
	const PolynomialVentura::Polynomial<8> t57 = B23*t37;
	const PolynomialVentura::Polynomial<8> t146 = B43*t23;
	const PolynomialVentura::Polynomial<8> t58 = t56 + t57 - t146;
	const PolynomialVentura::Polynomial<8> t59 = B53*t13;
	const PolynomialVentura::Polynomial<8> t60 = B33*t37;
	const PolynomialVentura::Polynomial<8> t102 = B43*t25;
	const PolynomialVentura::Polynomial<8> t61 = t59 + t60 - t102;
	const PolynomialVentura::Polynomial<12> t64 = B34*t35;
	const PolynomialVentura::Polynomial<12> t67 = B14*t55;
	const PolynomialVentura::Polynomial<12> t68 = t64 + t67 - B24*t45 - B54*t19;
	const PolynomialVentura::Polynomial<5> t69 = B11*B62;
	const PolynomialVentura::Polynomial<5> t76 = B12*B61;
	const PolynomialVentura::Polynomial<5> t70 = t69 - t76;
	const PolynomialVentura::Polynomial<5> t71 = B21*B62;
	const PolynomialVentura::Polynomial<5> t77 = B22*B61;
	const PolynomialVentura::Polynomial<5> t72 = t71 - t77;
	const PolynomialVentura::Polynomial<5> t73 = B31*B62;
	const PolynomialVentura::Polynomial<5> t83 = B32*B61;
	const PolynomialVentura::Polynomial<5> t74 = t73 - t83;
	const PolynomialVentura::Polynomial<8> t75 = B63*t3;
	const PolynomialVentura::Polynomial<8> t78 = B13*t72;
	const PolynomialVentura::Polynomial<8> t124 = B23*t70;
	const PolynomialVentura::Polynomial<8> t79 = t75 + t78 - t124;
	const PolynomialVentura::Polynomial<4> t80 = B51*B62;
	const PolynomialVentura::Polynomial<4> t87 = B52*B61;
	const PolynomialVentura::Polynomial<4> t81 = t80 - t87;
	const PolynomialVentura::Polynomial<8> t82 = B63*t5;
	const PolynomialVentura::Polynomial<8> t84 = B13*t74;
	const PolynomialVentura::Polynomial<8> t105 = B33*t70;
	const PolynomialVentura::Polynomial<8> t85 = t82 + t84 - t105;
	const PolynomialVentura::Polynomial<7> t86 = B63*t21;
	const PolynomialVentura::Polynomial<7> t88 = B13*t81;
	const PolynomialVentura::Polynomial<7> t109 = B53*t70;
	const PolynomialVentura::Polynomial<7> t89 = t86 + t88 - t109;
	const PolynomialVentura::Polynomial<8> t90 = B63*t9;
	const PolynomialVentura::Polynomial<8> t91 = B23*t74;
	const PolynomialVentura::Polynomial<8> t126 = B33*t72;
	const PolynomialVentura::Polynomial<8> t92 = t90 + t91 - t126;
	const PolynomialVentura::Polynomial<7> t93 = B63*t23;
	const PolynomialVentura::Polynomial<7> t94 = B23*t81;
	const PolynomialVentura::Polynomial<7> t150 = B53*t72;
	const PolynomialVentura::Polynomial<7> t95 = t93 + t94 - t150;
	const PolynomialVentura::Polynomial<7> t96 = B63*t25;
	const PolynomialVentura::Polynomial<7> t97 = B33*t81;
	const PolynomialVentura::Polynomial<7> t111 = B53*t74;
	const PolynomialVentura::Polynomial<7> t98 = t96 + t97 - t111;
	const PolynomialVentura::Polynomial<12> t100 = B44*t45;
	const PolynomialVentura::Polynomial<12> t103 = B14*t61;
	const PolynomialVentura::Polynomial<12> t104 = t100 + t103 - B34*t49 - B54*t41;
	const PolynomialVentura::Polynomial<5> t106 = B41*B62;
	const PolynomialVentura::Polynomial<5> t114 = B42*B61;
	const PolynomialVentura::Polynomial<5> t107 = t106 - t114;
	const PolynomialVentura::Polynomial<11> t108 = B64*t45;
	const PolynomialVentura::Polynomial<11> t110 = B34*t89;
	const PolynomialVentura::Polynomial<11> t112 = t108 + t110 - B14*t98 - B54*t85;
	const PolynomialVentura::Polynomial<8> t113 = B63*t7;
	const PolynomialVentura::Polynomial<8> t115 = B13*t107;
	const PolynomialVentura::Polynomial<8> t133 = B43*t70;
	const PolynomialVentura::Polynomial<8> t116 = t113 + t115 - t133;
	const PolynomialVentura::Polynomial<8> t117 = B63*t13;
	const PolynomialVentura::Polynomial<8> t118 = B33*t107;
	const PolynomialVentura::Polynomial<8> t136 = B43*t74;
	const PolynomialVentura::Polynomial<8> t119 = t117 + t118 - t136;
	const PolynomialVentura::Polynomial<7> t120 = B63*t37;
	const PolynomialVentura::Polynomial<7> t121 = B43*t81;
	const PolynomialVentura::Polynomial<7> t154 = B53*t107;
	const PolynomialVentura::Polynomial<7> t122 = t120 + t121 - t154;
	const PolynomialVentura::Polynomial<12> t123 = B64*t19;
	const PolynomialVentura::Polynomial<12> t125 = B24*t85;
	const PolynomialVentura::Polynomial<12> t127 = t123 + t125 - B14*t92 - B34*t79;
	const PolynomialVentura::Polynomial<13> t128 = B44*t19;
	const PolynomialVentura::Polynomial<13> t130 = B24*t41;
	const PolynomialVentura::Polynomial<13> t132 = t128 + t130 - B34*t30 - B14*t52;
	const PolynomialVentura::Polynomial<12> t134 = B64*t41;
	const PolynomialVentura::Polynomial<12> t135 = B34*t116;
	const PolynomialVentura::Polynomial<12> t137 = t134 + t135 - B44*t85 - B14*t119;
	const PolynomialVentura::Polynomial<8> t138 = B63*t11;
	const PolynomialVentura::Polynomial<8> t139 = B23*t107;
	const PolynomialVentura::Polynomial<8> t143 = B43*t72;
	const PolynomialVentura::Polynomial<8> t140 = t138 + t139 - t143;
	const PolynomialVentura::Polynomial<12> t141 = B64*t30;
	const PolynomialVentura::Polynomial<12> t142 = B24*t116;
	const PolynomialVentura::Polynomial<12> t144 = t141 + t142 - B44*t79 - B14*t140;
	const PolynomialVentura::Polynomial<12> t145 = B44*t35;
	const PolynomialVentura::Polynomial<12> t147 = B14*t58;
	const PolynomialVentura::Polynomial<11> t148 = B64*t35;
	const PolynomialVentura::Polynomial<11> t149 = B24*t89;
	const PolynomialVentura::Polynomial<11> t151 = t148 + t149 - B14*t95 - B54*t79;
	const PolynomialVentura::Polynomial<11> t152 = B64*t49;
	const PolynomialVentura::Polynomial<11> t153 = B44*t89;
	const PolynomialVentura::Polynomial<11> t155 = t152 + t153 - B14*t122 - B54*t116;
	const PolynomialVentura::Polynomial<12> t156 = B64*t52;
	const PolynomialVentura::Polynomial<12> t157 = B34*t140;
	const PolynomialVentura::Polynomial<12> t158 = t156 + t157 - B44*t92 - B24*t119;
	const PolynomialVentura::Polynomial<12> t159 = B44*t55;
	const PolynomialVentura::Polynomial<12> t160 = B24*t61;
	const PolynomialVentura::Polynomial<11> t161 = B64*t55;
	const PolynomialVentura::Polynomial<11> t162 = B34*t95;
	const PolynomialVentura::Polynomial<11> t163 = t161 + t162 - B24*t98 - B54*t92;
	const PolynomialVentura::Polynomial<11> t164 = B64*t58;
	const PolynomialVentura::Polynomial<11> t165 = B44*t95;
	const PolynomialVentura::Polynomial<11> t166 = t164 + t165 - B24*t122 - B54*t140;
	const PolynomialVentura::Polynomial<11> t167 = B64*t61;
	const PolynomialVentura::Polynomial<11> t168 = B44*t98;
	const PolynomialVentura::Polynomial<11> t169 = t167 + t168 - B34*t122 - B54*t119;
	const PolynomialVentura::Polynomial<20> detB = B66*(B45*t68 - B15*(t159 + t160 - B34*t58 - B54*t52) - B35*(t145 + t147 - B24*t49 - B54*t30) + B25*t104 + B55*t132) + B36*(B65*(t145 + t147 - B24*t49 - B54*t30) + B25*t155 - B15*t166 - B45*t151 + B55*t144) + B16*(B65*(t159 + t160 - B34*t58 - B54*t52) - B25*t169 + B35*t166 - B45*t163 + B55*t158) - B46*(B65*t68 + B25*t112 - B15*t163 + B55*t127 - B35*t151) - B26*(B65*t104 - B45*t112 - B15*t169 + B35*t155 + B55*t137) - B56*(B15*t158 - B45*t127 - B25*t137 + B35*t144 + B65*t132);

	std::vector<double> zsolns;
	detB.realRootsSturm(-15.*M_PI / 180., 15.*M_PI / 180., zsolns);

	rsolns.clear();
	rsolns.reserve(zsolns.size());
	for (size_t i = 0; i < zsolns.size(); i++)
	{
		Eigen::Matrix<double, 6, 6> Bz;
		Bz <<
			B11.eval(zsolns[i]), B12.eval(zsolns[i]), B13.eval(zsolns[i]), B14.eval(zsolns[i]), B15.eval(zsolns[i]), B16.eval(zsolns[i]),
			B21.eval(zsolns[i]), B22.eval(zsolns[i]), B23.eval(zsolns[i]), B24.eval(zsolns[i]), B25.eval(zsolns[i]), B26.eval(zsolns[i]),
			B31.eval(zsolns[i]), B32.eval(zsolns[i]), B33.eval(zsolns[i]), B34.eval(zsolns[i]), B35.eval(zsolns[i]), B36.eval(zsolns[i]),
			B41.eval(zsolns[i]), B42.eval(zsolns[i]), B43.eval(zsolns[i]), B44.eval(zsolns[i]), B45.eval(zsolns[i]), B46.eval(zsolns[i]),
			B51.eval(zsolns[i]), B52.eval(zsolns[i]), B53.eval(zsolns[i]), B54.eval(zsolns[i]), B55.eval(zsolns[i]), B56.eval(zsolns[i]),
			B61.eval(zsolns[i]), B62.eval(zsolns[i]), B63.eval(zsolns[i]), B64.eval(zsolns[i]), B65.eval(zsolns[i]), B66.eval(zsolns[i]);
		Eigen::Matrix<double, 6, 1> xysoln;
		xysoln << 0, 0, 0, 0, 0, 1;
		xysoln = Bz.householderQr().solve(xysoln);
		const double xsoln = xysoln(3) / xysoln(5);
		const double ysoln = xysoln(4) / xysoln(5);
		Eigen::Matrix<double, 3, 1> rsoln;
		rsoln << xsoln, ysoln, zsolns[i];
		rsolns.push_back(rsoln);
	}

}

namespace opengv
{
namespace relative_pose
{
namespace modules
{

struct Ge_step : OptimizationFunctor<double>
{
  const Eigen::Matrix3d & _xxF;
  const Eigen::Matrix3d & _yyF;
  const Eigen::Matrix3d & _zzF;
  const Eigen::Matrix3d & _xyF;
  const Eigen::Matrix3d & _yzF;
  const Eigen::Matrix3d & _zxF;
  const Eigen::Matrix<double,3,9> & _x1P;
  const Eigen::Matrix<double,3,9> & _y1P;
  const Eigen::Matrix<double,3,9> & _z1P;
  const Eigen::Matrix<double,3,9> & _x2P;
  const Eigen::Matrix<double,3,9> & _y2P;
  const Eigen::Matrix<double,3,9> & _z2P;
  const Eigen::Matrix<double,9,9> & _m11P;
  const Eigen::Matrix<double,9,9> & _m12P;
  const Eigen::Matrix<double,9,9> & _m22P;

  Ge_step(
    const Eigen::Matrix3d & xxF,
    const Eigen::Matrix3d & yyF,
    const Eigen::Matrix3d & zzF,
    const Eigen::Matrix3d & xyF,
    const Eigen::Matrix3d & yzF,
    const Eigen::Matrix3d & zxF,
    const Eigen::Matrix<double,3,9> & x1P,
    const Eigen::Matrix<double,3,9> & y1P,
    const Eigen::Matrix<double,3,9> & z1P,
    const Eigen::Matrix<double,3,9> & x2P,
    const Eigen::Matrix<double,3,9> & y2P,
    const Eigen::Matrix<double,3,9> & z2P,
    const Eigen::Matrix<double,9,9> & m11P,
    const Eigen::Matrix<double,9,9> & m12P,
    const Eigen::Matrix<double,9,9> & m22P ) :
    OptimizationFunctor<double>(3,3),
    _xxF(xxF),_yyF(yyF),_zzF(zzF),_xyF(xyF),_yzF(yzF),_zxF(zxF),
    _x1P(x1P),_y1P(y1P),_z1P(z1P),_x2P(x2P),_y2P(y2P),_z2P(z2P),
    _m11P(m11P),_m12P(m12P),_m22P(m22P) {}

  int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
  {
    cayley_t cayley = x;
    Eigen::Matrix<double,1,3> jacobian;
    ge::getCostWithJacobian(
        _xxF,_yyF,_zzF,_xyF,_yzF,_zxF,
        _x1P,_y1P,_z1P,_x2P,_y2P,_z2P,_m11P,_m12P,_m22P,cayley,jacobian,1);

    fvec[0] = jacobian(0,0);
    fvec[1] = jacobian(0,1);
    fvec[2] = jacobian(0,2);
    
    return 0;
  }
};

}
}
}

void
opengv::relative_pose::modules::ge_main(
    const Eigen::Matrix3d & xxF,
    const Eigen::Matrix3d & yyF,
    const Eigen::Matrix3d & zzF,
    const Eigen::Matrix3d & xyF,
    const Eigen::Matrix3d & yzF,
    const Eigen::Matrix3d & zxF,
    const Eigen::Matrix<double,3,9> & x1P,
    const Eigen::Matrix<double,3,9> & y1P,
    const Eigen::Matrix<double,3,9> & z1P,
    const Eigen::Matrix<double,3,9> & x2P,
    const Eigen::Matrix<double,3,9> & y2P,
    const Eigen::Matrix<double,3,9> & z2P,
    const Eigen::Matrix<double,9,9> & m11P,
    const Eigen::Matrix<double,9,9> & m12P,
    const Eigen::Matrix<double,9,9> & m22P,
    const cayley_t & startingPoint,
    geOutput_t & output )
{
  //this one doesn't work, probably because of double numerical differentiation
  //use ge_main2, which is an implementation of gradient descent
  const int n=3;
  VectorXd x(n);

  x = startingPoint;
  Ge_step functor(xxF,yyF,zzF,xyF,yzF,zxF,
      x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P);
  NumericalDiff<Ge_step> numDiff(functor);
  LevenbergMarquardt< NumericalDiff<Ge_step> > lm(numDiff);

  lm.resetParameters();
  lm.parameters.ftol = 0.000001;//1.E1*NumTraits<double>::epsilon();
  lm.parameters.xtol = 1.E1*NumTraits<double>::epsilon();
  lm.parameters.maxfev = 100;
  lm.minimize(x);

  cayley_t cayley = x;
  rotation_t R = math::cayley2rot(cayley);
  
  Eigen::Matrix4d G = ge::composeG(xxF,yyF,zzF,xyF,yzF,zxF,
      x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley);
  
  Eigen::EigenSolver< Eigen::Matrix4d > Eig(G,true);
  Eigen::Matrix<std::complex<double>,4,1> D_complex = Eig.eigenvalues();
  Eigen::Matrix<std::complex<double>,4,4> V_complex = Eig.eigenvectors();
  Eigen::Vector4d D;
  Eigen::Matrix4d V;
  for(size_t i = 0; i < 4; i++)
  {
    D[i] = D_complex[i].real();
    for(size_t j = 0; j < 4; j++)
      V(i,j) = V_complex(i,j).real();
  }

  double factor = V(3,0);
  Eigen::Vector4d t = (1.0/factor) * V.col(0);

  output.translation = t;
  output.rotation = R;
  output.eigenvalues = D;
  output.eigenvectors = V;
}

void
opengv::relative_pose::modules::ge_main2(
    const Eigen::Matrix3d & xxF,
    const Eigen::Matrix3d & yyF,
    const Eigen::Matrix3d & zzF,
    const Eigen::Matrix3d & xyF,
    const Eigen::Matrix3d & yzF,
    const Eigen::Matrix3d & zxF,
    const Eigen::Matrix<double,3,9> & x1P,
    const Eigen::Matrix<double,3,9> & y1P,
    const Eigen::Matrix<double,3,9> & z1P,
    const Eigen::Matrix<double,3,9> & x2P,
    const Eigen::Matrix<double,3,9> & y2P,
    const Eigen::Matrix<double,3,9> & z2P,
    const Eigen::Matrix<double,9,9> & m11P,
    const Eigen::Matrix<double,9,9> & m12P,
    const Eigen::Matrix<double,9,9> & m22P,
    const cayley_t & startingPoint,
    geOutput_t & output )
{
  //todo: the optimization strategy is something that can possibly be improved:
  //-one idea is to check the gradient at the new sampling point, if that derives
  // too much, we have to stop
  //-another idea consists of having linear change of lambda, instead of exponential (safer, but slower)
  
  double lambda = 0.01;
  double maxLambda = 0.08;
  double modifier = 2.0;
  int maxIterations = 50;
  double min_xtol = 0.00001;
  bool disablingIncrements = true;
  bool print = false;

  cayley_t cayley;
  
  double disturbanceAmplitude = 0.3;
  bool found = false;
  int randomTrialCount = 0;
  
  while( !found && randomTrialCount < 5 )
  {
    if(randomTrialCount > 2)
      disturbanceAmplitude = 0.6;
	
    if( randomTrialCount == 0 )
      cayley = startingPoint;
    else
    {
      cayley = startingPoint;
      Eigen::Vector3d disturbance;
      disturbance[0] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0*disturbanceAmplitude;
      disturbance[1] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0*disturbanceAmplitude;
      disturbance[2] = (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0*disturbanceAmplitude;
      cayley += disturbance;
    }
	
    lambda = 0.01;
    int iterations = 0;
    double smallestEV = ge::getCost(xxF,yyF,zzF,xyF,yzF,zxF,
        x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley,1);
    
    while( iterations < maxIterations )
    {
      Eigen::Matrix<double,1,3> jacobian;
      ge::getQuickJacobian(xxF,yyF,zzF,xyF,yzF,zxF,
          x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley,smallestEV,jacobian,1);
      
      double norm = sqrt(pow(jacobian[0],2.0) + pow(jacobian[1],2.0) + pow(jacobian[2],2.0));
      cayley_t normalizedJacobian = (1/norm) * jacobian.transpose();
      
      cayley_t samplingPoint = cayley - lambda * normalizedJacobian;
      double samplingEV = ge::getCost(xxF,yyF,zzF,xyF,yzF,zxF,
          x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,samplingPoint,1);
      
      if(print)
      {
        std::cout << iterations << ": " << samplingPoint.transpose();
        std::cout << " lambda: " << lambda << " EV: " << samplingEV << std::endl;
      }
      
      if( iterations == 0 || !disablingIncrements )
      {
        while( samplingEV < smallestEV )
        {
          smallestEV = samplingEV;
          if( lambda * modifier > maxLambda )
            break;
          lambda *= modifier;
          samplingPoint = cayley - lambda * normalizedJacobian;
          samplingEV = ge::getCost(xxF,yyF,zzF,xyF,yzF,zxF,
              x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,samplingPoint,1);
          
          if(print)
          {
            std::cout << iterations << ": " << samplingPoint.transpose();
            std::cout << " lambda: " << lambda << " EV: " << samplingEV << std::endl;
          }
        }
      }
      
      while( samplingEV > smallestEV )
      {
        lambda /= modifier;
        samplingPoint = cayley - lambda * normalizedJacobian;
        samplingEV = ge::getCost(xxF,yyF,zzF,xyF,yzF,zxF,
            x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,samplingPoint,1);
        
        if(print)
        {
          std::cout << iterations << ": " << samplingPoint.transpose();
          std::cout << " lambda: " << lambda << " EV: " << samplingEV << std::endl;
        }
      }
      
      //apply update
      cayley = samplingPoint;
      smallestEV = samplingEV;
      
      //stopping condition (check if the update was too small)
      if( lambda < min_xtol )
        break;
      
      iterations++;
    }
    
    //try to see if we can robustly identify each time we enter up in the wrong minimum
    if( cayley.norm() < 0.01 )
    {
      //we are close to the origin, test the EV 2
      double ev2 = ge::getCost(xxF,yyF,zzF,xyF,yzF,zxF,
            x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley,0);
      if( ev2 > 0.001 )
        randomTrialCount++;
      else
        found = true;
    }
    else
      found = true;
  }
  
  Eigen::Matrix4d G = ge::composeG(xxF,yyF,zzF,xyF,yzF,zxF,
      x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley);
  
  Eigen::EigenSolver< Eigen::Matrix4d > Eig(G,true);
  Eigen::Matrix<std::complex<double>,4,1> D_complex = Eig.eigenvalues();
  Eigen::Matrix<std::complex<double>,4,4> V_complex = Eig.eigenvectors();
  Eigen::Vector4d D;
  Eigen::Matrix4d V;
  for(size_t i = 0; i < 4; i++)
  {
    D[i] = D_complex[i].real();
    for(size_t j = 0; j < 4; j++)
      V(i,j) = V_complex(i,j).real();
  }

  double factor = V(3,0);
  Eigen::Vector4d t = (1.0/factor) * V.col(0);

  output.translation = t;
  output.rotation = math::cayley2rot(cayley);
  output.eigenvalues = D;
  output.eigenvectors = V;
}

void
opengv::relative_pose::modules::ge_plot(
    const Eigen::Matrix3d & xxF,
    const Eigen::Matrix3d & yyF,
    const Eigen::Matrix3d & zzF,
    const Eigen::Matrix3d & xyF,
    const Eigen::Matrix3d & yzF,
    const Eigen::Matrix3d & zxF,
    const Eigen::Matrix<double,3,9> & x1P,
    const Eigen::Matrix<double,3,9> & y1P,
    const Eigen::Matrix<double,3,9> & z1P,
    const Eigen::Matrix<double,3,9> & x2P,
    const Eigen::Matrix<double,3,9> & y2P,
    const Eigen::Matrix<double,3,9> & z2P,
    const Eigen::Matrix<double,9,9> & m11P,
    const Eigen::Matrix<double,9,9> & m12P,
    const Eigen::Matrix<double,9,9> & m22P,
    geOutput_t & output )
{
  cayley_t cayley = math::rot2cayley(output.rotation);
  
  ge::getEV(xxF,yyF,zzF,xyF,yzF,zxF,
      x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley,output.eigenvalues);
  
  output.eigenvectors = ge::composeG(xxF,yyF,zzF,xyF,yzF,zxF,
      x1P,y1P,z1P,x2P,y2P,z2P,m11P,m12P,m22P,cayley);
}
