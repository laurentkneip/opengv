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

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <opengv/math/Sturm.hpp>

#include "time_measurement.hpp"


using namespace std;
using namespace Eigen;
using namespace opengv;

int main( int argc, char** argv )
{
  std::vector<double> coeffs;
  coeffs.push_back(1.0);
  coeffs.push_back(4.0);
  coeffs.push_back(1.0);
  coeffs.push_back(-6.0);
  
  //timer
  struct timeval tic;
  struct timeval toc;
  size_t iterations = 50;
  
  //for now just construct the problem
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    math::Sturm sturm(coeffs);
  gettimeofday( &toc, 0 );
  double time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "the initialization takes " << time << " seconds" << std::endl;
  
  math::Sturm sturm(coeffs);
  
  //test the lagrangian bounds
  gettimeofday( &tic, 0 );
  double bound;
  for(size_t i = 0; i < iterations; i++)
    bound = sturm.computeLagrangianBound();
  gettimeofday( &toc, 0 );
  time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "the initial bound computation takes " << time << " seconds" << std::endl;
  std::cout << "the bound evaluates to " << bound << std::endl;
  
  //now evaluate the chain to verify that the number of roots is 3
  gettimeofday( &tic, 0 );
  size_t numberRoots;
  for(size_t i = 0; i < iterations; i++)
    numberRoots = sturm.evaluateChain(-bound) - sturm.evaluateChain(bound);
  gettimeofday( &toc, 0 );
  time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "the evaluation of two bounds takes " << time << " seconds" << std::endl;
  std::cout << "the obtained number of roots is " << numberRoots << std::endl;
  
  //now test the bracketing mechanism
  gettimeofday( &tic, 0 );
  std::vector<double> roots;
  for(size_t i = 0; i < iterations; i++)
  {
    roots.clear();
    sturm.bracketRoots(roots,0.5);
  }
  gettimeofday( &toc, 0 );
  time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "the bracketing of the roots took " << time << " seconds" << std::endl;
  std::cout << "the obtained brackets are:" << std::endl;
  std::vector<double>::iterator it = roots.begin();
  while( it != roots.end() )
  {
    std::cout << "root:" << (*it) << std::endl; 
    it++;
  }
  
  //now test the entire root-finding mechanism
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
    roots = sturm.findRoots();
  gettimeofday( &toc, 0 );
  time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "the entire finding of the roots took " << time << " seconds" << std::endl;
  std::cout << "the obtained roots are:" << std::endl;
  std::vector<double>::iterator it2 = roots.begin();
  while( it2 != roots.end() )
  {
    std::cout << (*it2) << std::endl;
    it2++;
  }
  
  //now test the new root-finding with inbuild gradient descent
  gettimeofday( &tic, 0 );
  for(size_t i = 0; i < iterations; i++)
  {
    roots.clear();
    sturm.findRoots2(roots);
  }
  gettimeofday( &toc, 0 );
  time = TIMETODOUBLE(timeval_minus(toc,tic)) / iterations;
  
  std::cout << "the entire finding of the roots took " << time << " seconds" << std::endl;
  std::cout << "the obtained roots are:" << std::endl;
  it2 = roots.begin();
  while( it2 != roots.end() )
  {
    std::cout << (*it2) << std::endl;
    it2++;
  }
  
  return 0;
}
