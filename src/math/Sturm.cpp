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


#include <opengv/math/Sturm.hpp>
#include <iostream>

opengv::math::Sturm::Sturm( const Eigen::MatrixXd & p ) :
    _C(Eigen::MatrixXd(p.cols(),p.cols()))
{
  _dimension = (size_t) _C.cols();
  _C = Eigen::MatrixXd::Zero(_dimension,_dimension);
  _C.row(0) = p;

  for( size_t i = 1; i < _dimension; i++ )
    _C(1,i) = _C(0,i-1) * (_dimension-i);

  for( size_t i = 2; i < _dimension; i++ )
  {
    Eigen::MatrixXd p1(1,_dimension-(i-2));
    p1 = _C.block(i-2,i-2,1,_dimension-(i-2));
    Eigen::MatrixXd p2(1,_dimension-(i-1));
    p2 = _C.block(i-1,i-1,1,_dimension-(i-1));
    Eigen::MatrixXd r(1,_dimension-i);
    computeNegatedRemainder(p1,p2,r);
    _C.block(i,i,1,_dimension-i) = r;
  }
}

opengv::math::Sturm::Sturm( const std::vector<double> & p ) :
    _C(Eigen::MatrixXd(p.size(),p.size()))
{
  _dimension = (size_t) _C.cols();
  _C = Eigen::MatrixXd::Zero(_dimension,_dimension);

  for( size_t i = 0; i < _dimension; i++ )
    _C(0,i) = p[i];

  for( size_t i = 1; i < _dimension; i++ )
    _C(1,i) = _C(0,i-1) * (_dimension-i);

  for( size_t i = 2; i < _dimension; i++ )
  {
    Eigen::MatrixXd p1(1,_dimension-(i-2));
    p1 = _C.block(i-2,i-2,1,_dimension-(i-2));
    Eigen::MatrixXd p2(1,_dimension-(i-1));
    p2 = _C.block(i-1,i-1,1,_dimension-(i-1));
    Eigen::MatrixXd r(1,_dimension-i);
    computeNegatedRemainder(p1,p2,r);
    _C.block(i,i,1,_dimension-i) = r;
  }
}

opengv::math::Sturm::~Sturm() {};

std::vector<double>
opengv::math::Sturm::findRoots()
{
  //bracket the roots
  std::vector<std::pair<double,double> > brackets = bracketRoots();

  //pollish
  std::vector<double> roots;
  roots.resize(brackets.size());
  Eigen::MatrixXd monomials(_dimension,1);
  monomials(_dimension-1,0) = 1.0;

  for(size_t r = 0; r < brackets.size(); r++ )
  {
    roots[r] = 0.5 * (brackets[r].first + brackets[r].second);

    //Now do Gauss iterations here
    //evaluate all monomials at the left bound
    for(size_t k = 0; k < 100; k++ )
    {
      for(size_t i = 2; i <= _dimension; i++)
        monomials(_dimension-i,0) = monomials(_dimension-i+1,0)*roots[r];

      Eigen::MatrixXd value = _C.row(0) * monomials;
      Eigen::MatrixXd derivative = _C.row(1) * monomials;

      //correction
      roots[r] = roots[r] - (value(0,0)/derivative(0,0));
    }
  }

  return roots;
}

std::vector<std::pair<double,double> >
opengv::math::Sturm::bracketRoots()
{
  double absoluteRightBound = computeLagrangianBound();
  double absoluteLeftBound = -absoluteRightBound;
  size_t totalNumberOfRoots =
      evaluateChain(absoluteLeftBound,absoluteRightBound);

  std::vector< std::pair<double,double> > brackets;
  brackets.reserve(totalNumberOfRoots);
  brackets.push_back(
      std::pair<double,double>(absoluteLeftBound,absoluteRightBound));
  std::vector<size_t> roots;
  roots.push_back(totalNumberOfRoots);

  while(brackets.size() < totalNumberOfRoots)
  {
    size_t numberBrackets = brackets.size();

    //Bisectioning of each bracket that needs to
    for(size_t i = 0; i < numberBrackets; i++)
    {
      //derive the bound of the left side
      double leftBound = brackets[i].first;
      double rightBound = 0.5*(leftBound + brackets[i].second);

      size_t numberRoots = evaluateChain(leftBound,rightBound);

      if(numberRoots == roots[i])
      {
        //the roots are all in the left part, just change the right bracket
        brackets[i].second = rightBound;
      }
      else
      {
        if(numberRoots == 0)
        {
          //the roots are all in the right part, just change the left bracket
          brackets[i].first = rightBound;
        }
        else
        {
          //split up the bracket
          brackets.push_back(
              std::pair<double,double>(rightBound,brackets[i].second));
          brackets[i].second = rightBound;
          roots.push_back(roots[i]-numberRoots);
          roots[i] = numberRoots;
        }
      }
    }
  }

  //Now all the roots should be isolated, we do another k steps per root to
  //tighten the brackets
  for(size_t i=0; i < totalNumberOfRoots; i++)
  {
    for(size_t k = 0; k < 10; k++)
    {
      //derive the bound of the left side
      double leftBound = brackets[i].first;
      double rightBound = 0.5*(leftBound + brackets[i].second);
      size_t numberRoots = evaluateChain(leftBound,rightBound);

      if(numberRoots == 0)
        brackets[i].first = rightBound;
      else
        brackets[i].second = rightBound;
    }
  }

  return brackets;
}

size_t
opengv::math::Sturm::evaluateChain(
    double leftBound,
    double rightBound )
{
  Eigen::MatrixXd monomials(_dimension,2);
  monomials(_dimension-1,0) = 1.0;
  monomials(_dimension-1,1) = 1.0;

  //evaluate all monomials at the bounds
  for(size_t i = 2; i <= _dimension; i++)
  {
    monomials(_dimension-i,0) = monomials(_dimension-i+1,0)*leftBound;
    monomials(_dimension-i,1) = monomials(_dimension-i+1,1)*rightBound;
  }

  //evaluate the sign of the polynomials at the bounds
  Eigen::MatrixXd signs(_dimension,2);
  signs = _C*monomials;

  //count the sign changes
  for(size_t i = 0; i < _dimension; i++)
  {
    signs(i,0) = signs(i,0) / fabs(signs(i,0));
    signs(i,1) = signs(i,1) / fabs(signs(i,1));
  }

  Eigen::MatrixXd temp(_dimension-1,2);
  temp = signs.block(0,0,_dimension-1,2) + signs.block(1,0,_dimension-1,2);
  size_t leftSignChanges = 0;
  size_t rightSignChanges = 0;
  for(size_t i= 0; i < _dimension-1; i++)
  {
    if( fabs(temp(i,0)) < 1.0 )
      leftSignChanges++;
    if( fabs(temp(i,1)) < 1.0 )
      rightSignChanges++;
  }

  size_t numberOfRoots = leftSignChanges - rightSignChanges;
  return numberOfRoots;
}

void
opengv::math::Sturm::computeNegatedRemainder(
    const Eigen::MatrixXd & p1,
    const Eigen::MatrixXd & p2,
    Eigen::MatrixXd & r )
{
  //we have to create 3 subtraction polynomials
  Eigen::MatrixXd poly_1(1,p1.cols());
  poly_1 = Eigen::MatrixXd::Zero(1,p1.cols());
  poly_1.block(0,0,1,p2.cols()) = (p1(0,0)/p2(0,0)) * p2;

  Eigen::MatrixXd poly_2(1,p1.cols());
  poly_2 = Eigen::MatrixXd::Zero(1,p1.cols());
  poly_2.block(0,1,1,p2.cols()) = (p1(0,1)/p2(0,0)) * p2;

  Eigen::MatrixXd poly_3(1,p1.cols());
  poly_3 = Eigen::MatrixXd::Zero(1,p1.cols());
  poly_3.block(0,1,1,p2.cols()) = (-p2(0,1)*p1(0,0)/pow(p2(0,0),2)) * p2;

  //compute remainder
  Eigen::MatrixXd remainder(1,p1.cols());
  remainder = p1 - poly_1 - poly_2 - poly_3;
  r = -remainder.block(0,2,1,p1.cols()-2);
}

double
opengv::math::Sturm::computeLagrangianBound()
{
  std::vector<double> coefficients;
  coefficients.reserve(_dimension-1);

  for(size_t i=0; i < _dimension-1; i++)
    coefficients.push_back(pow(fabs(_C(0,i+1)/_C(0,0)),(1.0/(i+1))));

  size_t j = 0;
  double max1 = -1.0;
  for( size_t i = 0; i < coefficients.size(); i++ )
  {
    if( coefficients[i] > max1 )
    {
      j = i;
      max1 = coefficients[i];
    }
  }

  coefficients[j] = -1.0;

  double max2 = -1.0;
  for( size_t i=0; i < coefficients.size(); i++ )
  {
    if( coefficients[i] > max2 )
      max2 = coefficients[i];
  }

  double bound = max1 + max2;
  return bound;
}
