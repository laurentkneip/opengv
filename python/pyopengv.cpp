#include <boost/python.hpp>
#include <iostream>
#include <vector>
#include <opengv/absolute_pose/AbsoluteAdapterBase.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include "types.hpp"

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#if (PY_VERSION_HEX < 0x03000000)
static void numpy_import_array_wrapper()
#else
static int* numpy_import_array_wrapper()
#endif
{
  /* Initialise numpy API and use 2/3 compatible return */
  import_array();
}




namespace pyopengv {

namespace bp = boost::python;
namespace bpn = boost::python::numeric;

typedef PyArrayContiguousView<double> pyarray_t;


opengv::bearingVector_t bearingVectorFromArray(
    const pyarray_t &array,
    size_t index )
{
  opengv::bearingVector_t v;
  v[0] = array.get(index, 0);
  v[1] = array.get(index, 1);
  v[2] = array.get(index, 2);
    for (int i = 0; i < 3; ++i) {
      std::cout << v[i] << " ";
    }
    std::cout << "  bearing\n";
  return v;
}


opengv::point_t pointFromArray(
    const pyarray_t &array,
    size_t index )
{
  opengv::point_t p;
  p[0] = array.get(index, 0);
  p[1] = array.get(index, 1);
  p[2] = array.get(index, 2);
    for (int i = 0; i < 3; ++i) {
      std::cout << p[i] << " ";
    }
    std::cout << "  point\n";
  return p;
}


bp::list listFromTransformations( const opengv::transformations_t &t )
{
  bp::list retn;
  for (size_t i = 0; i < t.size(); ++i) {
    std::cout << t[i] << "\n";
    Eigen::Matrix<double, 3, 4, Eigen::RowMajor> t_row_major = t[i];
    npy_intp shape[2] = {3, 4};
    retn.append(bpn_array_from_data(2, shape, t_row_major.data()));
  }
  return retn;
}


class CentralAbsoluteAdapter : public opengv::absolute_pose::AbsoluteAdapterBase
{
protected:
  using AbsoluteAdapterBase::_t;
  using AbsoluteAdapterBase::_R;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  /**
   * \brief Constructor. See protected class-members to understand parameters
   */
  CentralAbsoluteAdapter(
      bpn::array & bearingVectors,
      bpn::array & points )
    : _bearingVectors(bearingVectors)
    , _points(points)
  {}

  /**
   * \brief Constructor. See protected class-members to understand parameters
   */
  CentralAbsoluteAdapter(
      bpn::array & bearingVectors,
      bpn::array & points,
      bpn::array & R )
    : _bearingVectors(bearingVectors)
    , _points(points)
  {
    pyarray_t R_view(R);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        _R(i, j) = R_view.get(i, j);
      }
    }
  }

  /**
   * \brief Constructor. See protected class-members to understand parameters
   */
  CentralAbsoluteAdapter(
      bpn::array & bearingVectors,
      bpn::array & points,
      bpn::array & t,
      bpn::array & R )
    : _bearingVectors(bearingVectors)
    , _points(points)
  {
    pyarray_t t_view(t);
    for (int i = 0; i < 3; ++i) {
      _t(i) = t_view.get(i);
    }
    pyarray_t R_view(R);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        _R(i, j) = R_view.get(i, j);
      }
    }
  }
  /**
   * Destructor
   */
  virtual ~CentralAbsoluteAdapter() {}

  //Access of correspondences

  /** See parent-class */
  virtual opengv::bearingVector_t getBearingVector( size_t index ) const {
    opengv::bearingVector_t b = bearingVectorFromArray(_bearingVectors, index);
    return b;
  }
  /** See parent-class */
  virtual double getWeight( size_t index ) const {
    return 1.0;
  }
  /** See parent-class. Returns zero for this adapter. */
  virtual opengv::translation_t getCamOffset( size_t index ) const {
    return Eigen::Vector3d::Zero();
  }
  /** See parent-class Returns identity for this adapter. */
  virtual opengv::rotation_t getCamRotation( size_t index ) const {
    return opengv::rotation_t::Identity();
  }
  /** See parent-class */
  virtual opengv::point_t getPoint( size_t index ) const {
    return pointFromArray(_points, index);
  }
  /** See parent-class */
  virtual size_t getNumberCorrespondences() const {
    return _bearingVectors.shape(0);
  }

protected:
  /** Reference to the bearing-vectors expressed in the camera-frame */
  pyarray_t _bearingVectors;
  /** Reference to the points expressed in the world-frame. */
  pyarray_t _points;
};



bp::object absolute_pose_p3p_kneip( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter( v, p );

  opengv::transformations_t transformations =
      opengv::absolute_pose::p3p_kneip( adapter, 0, 1, 2 );

  return listFromTransformations(transformations);
}

}

BOOST_PYTHON_MODULE(pyopengv) {
  using namespace boost::python;

  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  numpy_import_array_wrapper();

  def("absolute_pose_p3p_kneip", pyopengv::absolute_pose_p3p_kneip);

}
