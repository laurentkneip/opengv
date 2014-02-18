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

/**
 * \file Indices.hpp
 * \brief Generic functor base for use with the Eigen-nonlinear optimization
 *        toolbox.
 */

#ifndef OPENGV_INDICES_HPP_
#define OPENGV_INDICES_HPP_

#include <stdlib.h>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

/**
 * \brief The namespace of this library.
 */
namespace opengv
{

/**
 * Index class used internally so we can efficiently access either all or
 * just a subset of the indices, using the same syntax
 */
struct Indices
{
  bool _useIndices;
  const std::vector<int> * _indices;
  size_t _numberCorrespondences;

  Indices( const std::vector<int> & indices) :
    _useIndices(true),
    _indices(&indices),
    _numberCorrespondences(indices.size())
  {}

  Indices(size_t numberCorrespondences) :
      _useIndices(false),
      _numberCorrespondences(numberCorrespondences)
  {}

  size_t size() const
  {
    return _numberCorrespondences;
  }

  int operator[](int i) const
  {
    if( _useIndices )
      return (*_indices)[i];
    return i;
  }
};

}

#endif /* OPENGV_INDICES_HPP_ */
