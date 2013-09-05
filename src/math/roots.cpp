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


#include <opengv/math/roots.hpp>
#include <complex>

std::vector<double>
opengv::math::o4_roots( const Eigen::MatrixXd & p )
{
  double A = p(0,0);
  double B = p(1,0);
  double C = p(2,0);
  double D = p(3,0);
  double E = p(4,0);

  double A_pw2 = A*A;
  double B_pw2 = B*B;
  double A_pw3 = A_pw2*A;
  double B_pw3 = B_pw2*B;
  double A_pw4 = A_pw3*A;
  double B_pw4 = B_pw3*B;

  double alpha = -3*B_pw2/(8*A_pw2)+C/A;
  double beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
  double gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

  double alpha_pw2 = alpha*alpha;
  double alpha_pw3 = alpha_pw2*alpha;

  std::complex<double> P (-alpha_pw2/12-gamma,0);
  std::complex<double> Q (-alpha_pw3/108+alpha*gamma/3-pow(beta,2)/8,0);
  std::complex<double> R = -Q/2.0+sqrt(pow(Q,2.0)/4.0+pow(P,3.0)/27.0);

  std::complex<double> U = pow(R,(1.0/3.0));
  std::complex<double> y;

  if (U.real() == 0)
    y = -5.0*alpha/6.0-pow(Q,(1.0/3.0));
  else
    y = -5.0*alpha/6.0-P/(3.0*U)+U;

  std::complex<double> w = sqrt(alpha+2.0*y);

  std::vector<double> realRoots;
  std::complex<double> temp;
  temp = -B/(4.0*A) + 0.5*(w+sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(w-sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(-w+sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(-w-sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  realRoots.push_back(temp.real());

  return realRoots;
}

std::vector<double>
opengv::math::o4_roots( const std::vector<double> & p )
{
  double A = p[0];
  double B = p[1];
  double C = p[2];
  double D = p[3];
  double E = p[4];

  double A_pw2 = A*A;
  double B_pw2 = B*B;
  double A_pw3 = A_pw2*A;
  double B_pw3 = B_pw2*B;
  double A_pw4 = A_pw3*A;
  double B_pw4 = B_pw3*B;

  double alpha = -3*B_pw2/(8*A_pw2)+C/A;
  double beta = B_pw3/(8*A_pw3)-B*C/(2*A_pw2)+D/A;
  double gamma = -3*B_pw4/(256*A_pw4)+B_pw2*C/(16*A_pw3)-B*D/(4*A_pw2)+E/A;

  double alpha_pw2 = alpha*alpha;
  double alpha_pw3 = alpha_pw2*alpha;

  std::complex<double> P (-alpha_pw2/12-gamma,0);
  std::complex<double> Q (-alpha_pw3/108+alpha*gamma/3-pow(beta,2)/8,0);
  std::complex<double> R = -Q/2.0+sqrt(pow(Q,2.0)/4.0+pow(P,3.0)/27.0);

  std::complex<double> U = pow(R,(1.0/3.0));
  std::complex<double> y;

  if (U.real() == 0)
    y = -5.0*alpha/6.0-pow(Q,(1.0/3.0));
  else
    y = -5.0*alpha/6.0-P/(3.0*U)+U;

  std::complex<double> w = sqrt(alpha+2.0*y);

  std::vector<double> realRoots;
  std::complex<double> temp;
  temp = -B/(4.0*A) + 0.5*(w+sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(w-sqrt(-(3.0*alpha+2.0*y+2.0*beta/w)));
  realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(-w+sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  realRoots.push_back(temp.real());
  temp = -B/(4.0*A) + 0.5*(-w-sqrt(-(3.0*alpha+2.0*y-2.0*beta/w)));
  realRoots.push_back(temp.real());

  return realRoots;
}
