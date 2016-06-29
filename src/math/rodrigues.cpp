/******************************************************************************
* Author:   Steffen Urban                                                    *
* Contact:  urbste@googlemail.com                                            *
* License:  Copyright (c) 2016 Steffen Urban, ANU. All rights reserved.      *
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

#include <opengv/math/rodrigues.hpp>

opengv::rotation_t
opengv::math::rodrigues2rot( const cayley_t & omega)
{
  rotation_t R = Eigen::Matrix3d::Identity();
	
  Eigen::Matrix3d skewW;
  skewW<<0.0, -omega(2), omega(1),
	  omega(2), 0.0, -omega(0),
	  -omega(1), omega(0), 0.0;

  double omega_norm = omega.norm();

  if (omega_norm > std::numeric_limits<double>::epsilon())
	  R = R + sin(omega_norm) / omega_norm*skewW 
	  + (1 - cos(omega_norm)) / (omega_norm*omega_norm)*(skewW*skewW);

  return R;
}


opengv::rodrigues_t
opengv::math::rot2rodrigues( const rotation_t & R )
{
  rodrigues_t omega;
  omega << 0.0, 0.0, 0.0;
  
  double trace = R.trace() - 1.0;
  double wnorm = acos(trace / 2.0);
  if (wnorm > std::numeric_limits<double>::epsilon())
  {
	  omega[0] = (R(2, 1) - R(1, 2));
	  omega[1] = (R(0, 2) - R(2, 0));
	  omega[2] = (R(1, 0) - R(0, 1));
	  double sc = wnorm / (2.0*sin(wnorm));
	  omega *= sc;
  }
  return omega;
}
