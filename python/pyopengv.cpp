#include <boost/python.hpp>
#include <iostream>
#include <vector>
#include <opengv/absolute_pose/AbsoluteAdapterBase.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/relative_pose/RelativeAdapterBase.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp>
#include <opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp>
#include <opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp>
#include <opengv/triangulation/methods.hpp>

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
  return p;
}

bp::object arrayFromPoints( const opengv::points_t &points )
{
  std::vector<double> data(points.size() * 3);
  for (size_t i = 0; i < points.size(); ++i) {
    data[3 * i + 0] = points[i][0];
    data[3 * i + 1] = points[i][1];
    data[3 * i + 2] = points[i][2];
  }
  npy_intp shape[2] = {(npy_intp)points.size(), 3};
  return bpn_array_from_data(2, shape, &data[0]);
}

bp::object arrayFromTranslation( const opengv::translation_t &t )
{
  npy_intp shape[1] = {3};
  return bpn_array_from_data(1, shape, t.data());
}

bp::object arrayFromRotation( const opengv::rotation_t &R )
{
  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> R_row_major = R;
  npy_intp shape[2] = {3, 3};
  return bpn_array_from_data(2, shape, R_row_major.data());
}

bp::list listFromRotations( const opengv::rotations_t &Rs )
{
  bp::list retn;
  for (size_t i = 0; i < Rs.size(); ++i) {
    retn.append(arrayFromRotation(Rs[i]));
  }
  return retn;
}

bp::object arrayFromEssential( const opengv::essential_t &E )
{
  Eigen::Matrix<double, 3, 3, Eigen::RowMajor> E_row_major = E;
  npy_intp shape[2] = {3, 3};
  return bpn_array_from_data(2, shape, E_row_major.data());
}

bp::list listFromEssentials( const opengv::essentials_t &Es )
{
  bp::list retn;
  for (size_t i = 0; i < Es.size(); ++i) {
    retn.append(arrayFromEssential(Es[i]));
  }
  return retn;
}

bp::object arrayFromTransformation( const opengv::transformation_t &t )
{
  Eigen::Matrix<double, 3, 4, Eigen::RowMajor> t_row_major = t;
  npy_intp shape[2] = {3, 4};
  return bpn_array_from_data(2, shape, t_row_major.data());
}

bp::list listFromTransformations( const opengv::transformations_t &t )
{
  bp::list retn;
  for (size_t i = 0; i < t.size(); ++i) {
    retn.append(arrayFromTransformation(t[i]));
  }
  return retn;
}

std::vector<int> getNindices( int n )
{
  std::vector<int> indices;
  for(int i = 0; i < n; i++)
    indices.push_back(i);
  return indices;
}


namespace absolute_pose {

class CentralAbsoluteAdapter : public opengv::absolute_pose::AbsoluteAdapterBase
{
protected:
  using AbsoluteAdapterBase::_t;
  using AbsoluteAdapterBase::_R;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  CentralAbsoluteAdapter(
      bpn::array & bearingVectors,
      bpn::array & points )
    : _bearingVectors(bearingVectors)
    , _points(points)
  {}

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

  virtual ~CentralAbsoluteAdapter() {}

  //Access of correspondences

  virtual opengv::bearingVector_t getBearingVector( size_t index ) const {
    return bearingVectorFromArray(_bearingVectors, index);
  }

  virtual double getWeight( size_t index ) const {
    return 1.0;
  }

  virtual opengv::translation_t getCamOffset( size_t index ) const {
    return Eigen::Vector3d::Zero();
  }

  virtual opengv::rotation_t getCamRotation( size_t index ) const {
    return opengv::rotation_t::Identity();
  }

  virtual opengv::point_t getPoint( size_t index ) const {
    return pointFromArray(_points, index);
  }

  virtual size_t getNumberCorrespondences() const {
    return _bearingVectors.shape(0);
  }

protected:
  pyarray_t _bearingVectors;
  pyarray_t _points;
};



bp::object p2p( bpn::array &v, bpn::array &p, bpn::array &R )
{
  CentralAbsoluteAdapter adapter(v, p, R);
  return arrayFromTranslation(
    opengv::absolute_pose::p2p(adapter, 0, 1));
}

bp::object p3p_kneip( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter(v, p);
  return listFromTransformations(
    opengv::absolute_pose::p3p_kneip(adapter, 0, 1, 2));
}

bp::object p3p_gao( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter(v, p);
  return listFromTransformations(
    opengv::absolute_pose::p3p_gao(adapter, 0, 1, 2));
}

bp::object gp3p( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter(v, p);
  return listFromTransformations(
    opengv::absolute_pose::gp3p(adapter, 0, 1, 2));
}

bp::object epnp( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter(v, p);
  return arrayFromTransformation(
    opengv::absolute_pose::epnp(adapter));
}

bp::object gpnp( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter(v, p);
  return arrayFromTransformation(
    opengv::absolute_pose::gpnp(adapter));
}

bp::object upnp( bpn::array &v, bpn::array &p )
{
  CentralAbsoluteAdapter adapter(v, p);
  return listFromTransformations(
    opengv::absolute_pose::upnp(adapter));
}

bp::object optimize_nonlinear( bpn::array &v,
                               bpn::array &p,
                               bpn::array &t,
                               bpn::array &R )
{
  CentralAbsoluteAdapter adapter(v, p, t, R);
  return arrayFromTransformation(
    opengv::absolute_pose::optimize_nonlinear(adapter));
}

bp::object ransac(
    bpn::array &v,
    bpn::array &p,
    std::string algo_name,
    double threshold,
    int max_iterations )
{
  using namespace opengv::sac_problems::absolute_pose;

  CentralAbsoluteAdapter adapter(v, p);

  // Create a ransac problem
  AbsolutePoseSacProblem::algorithm_t algorithm = AbsolutePoseSacProblem::KNEIP;
  if (algo_name == "TWOPT") algorithm = AbsolutePoseSacProblem::TWOPT;
  else if (algo_name == "KNEIP") algorithm = AbsolutePoseSacProblem::KNEIP;
  else if (algo_name == "GAO") algorithm = AbsolutePoseSacProblem::GAO;
  else if (algo_name == "EPNP") algorithm = AbsolutePoseSacProblem::EPNP;
  else if (algo_name == "GP3P") algorithm = AbsolutePoseSacProblem::GP3P;

  std::shared_ptr<AbsolutePoseSacProblem>
      absposeproblem_ptr(
        new AbsolutePoseSacProblem(adapter, algorithm));

  // Create a ransac solver for the problem
  opengv::sac::Ransac<AbsolutePoseSacProblem> ransac;

  ransac.sac_model_ = absposeproblem_ptr;
  ransac.threshold_ = threshold;
  ransac.max_iterations_ = max_iterations;

  // Solve
  ransac.computeModel();
  return arrayFromTransformation(ransac.model_coefficients_);
}



} // namespace absolute_pose


namespace relative_pose
{

class CentralRelativeAdapter : public opengv::relative_pose::RelativeAdapterBase
{
protected:
  using RelativeAdapterBase::_t12;
  using RelativeAdapterBase::_R12;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  CentralRelativeAdapter(
      bpn::array & bearingVectors1,
      bpn::array & bearingVectors2 )
    : _bearingVectors1(bearingVectors1)
    , _bearingVectors2(bearingVectors2)
  {}

  CentralRelativeAdapter(
      bpn::array & bearingVectors1,
      bpn::array & bearingVectors2,
      bpn::array & R12 )
    : _bearingVectors1(bearingVectors1)
    , _bearingVectors2(bearingVectors2)
  {
    pyarray_t R12_view(R12);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        _R12(i, j) = R12_view.get(i, j);
      }
    }
  }

  CentralRelativeAdapter(
      bpn::array & bearingVectors1,
      bpn::array & bearingVectors2,
      bpn::array & t12,
      bpn::array & R12 )
    : _bearingVectors1(bearingVectors1)
    , _bearingVectors2(bearingVectors2)
  {
    pyarray_t t12_view(t12);
    for (int i = 0; i < 3; ++i) {
      _t12(i) = t12_view.get(i);
    }
    pyarray_t R12_view(R12);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        _R12(i, j) = R12_view.get(i, j);
      }
    }
  }

  virtual ~CentralRelativeAdapter() {}

  virtual opengv::bearingVector_t getBearingVector1( size_t index ) const {
    return bearingVectorFromArray(_bearingVectors1, index);
  }

  virtual opengv::bearingVector_t getBearingVector2( size_t index ) const {
    return bearingVectorFromArray(_bearingVectors2, index);
  }

  virtual double getWeight( size_t index ) const {
    return 1.0;
  }

  virtual opengv::translation_t getCamOffset1( size_t index ) const {
    return Eigen::Vector3d::Zero();
  }

  virtual opengv::rotation_t getCamRotation1( size_t index ) const {
    return opengv::rotation_t::Identity();
  }

  virtual opengv::translation_t getCamOffset2( size_t index ) const {
    return Eigen::Vector3d::Zero();
  }

  virtual opengv::rotation_t getCamRotation2( size_t index ) const {
    return opengv::rotation_t::Identity();
  }

  virtual size_t getNumberCorrespondences() const {
    return _bearingVectors1.shape(0);
  }

protected:
  pyarray_t _bearingVectors1;
  pyarray_t _bearingVectors2;
};


bp::object twopt( bpn::array &b1, bpn::array &b2, bpn::array &R )
{
  CentralRelativeAdapter adapter(b1, b2, R);
  return arrayFromTranslation(
    opengv::relative_pose::twopt(adapter, true, 0, 1));
}

bp::object twopt_rotationOnly( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return arrayFromRotation(
    opengv::relative_pose::twopt_rotationOnly(adapter, 0, 1));
}

bp::object rotationOnly( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return arrayFromRotation(
    opengv::relative_pose::rotationOnly(adapter));
}

bp::object fivept_nister( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return listFromEssentials(
    opengv::relative_pose::fivept_nister(adapter));
}

bp::object fivept_kneip( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return listFromRotations(
    opengv::relative_pose::fivept_kneip(adapter, getNindices(5)));
}

bp::object sevenpt( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return listFromEssentials(
    opengv::relative_pose::sevenpt(adapter));
}

bp::object eightpt( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return arrayFromEssential(
    opengv::relative_pose::eightpt(adapter));
}

bp::object eigensolver( bpn::array &b1, bpn::array &b2, bpn::array &R )
{
  CentralRelativeAdapter adapter(b1, b2, R);
  return arrayFromRotation(
    opengv::relative_pose::eigensolver(adapter));
}

bp::object sixpt( bpn::array &b1, bpn::array &b2 )
{
  CentralRelativeAdapter adapter(b1, b2);
  return listFromRotations(
    opengv::relative_pose::sixpt(adapter));
}

bp::object optimize_nonlinear( bpn::array &b1,
                               bpn::array &b2,
                               bpn::array &t12,
                               bpn::array &R12 )
{
  CentralRelativeAdapter adapter(b1, b2, t12, R12);
  return arrayFromTransformation(
    opengv::relative_pose::optimize_nonlinear(adapter));
}

bp::object ransac(
    bpn::array &b1,
    bpn::array &b2,
    std::string algo_name,
    double threshold,
    int max_iterations )
{
  using namespace opengv::sac_problems::relative_pose;

  CentralRelativeAdapter adapter(b1, b2);

  // Create a ransac problem
  CentralRelativePoseSacProblem::algorithm_t algorithm = CentralRelativePoseSacProblem::NISTER;
  if (algo_name == "STEWENIUS") algorithm = CentralRelativePoseSacProblem::STEWENIUS;
  else if (algo_name == "NISTER") algorithm = CentralRelativePoseSacProblem::NISTER;
  else if (algo_name == "SEVENPT") algorithm = CentralRelativePoseSacProblem::SEVENPT;
  else if (algo_name == "EIGHTPT") algorithm = CentralRelativePoseSacProblem::EIGHTPT;

  std::shared_ptr<CentralRelativePoseSacProblem>
      relposeproblem_ptr(
        new CentralRelativePoseSacProblem(adapter, algorithm));

  // Create a ransac solver for the problem
  opengv::sac::Ransac<CentralRelativePoseSacProblem> ransac;

  ransac.sac_model_ = relposeproblem_ptr;
  ransac.threshold_ = threshold;
  ransac.max_iterations_ = max_iterations;

  // Solve
  ransac.computeModel();
  return arrayFromTransformation(ransac.model_coefficients_);
}

bp::object ransac_rotationOnly(
    bpn::array &b1,
    bpn::array &b2,
    double threshold,
    int max_iterations )
{
  using namespace opengv::sac_problems::relative_pose;

  CentralRelativeAdapter adapter(b1, b2);

  std::shared_ptr<RotationOnlySacProblem>
      relposeproblem_ptr(
        new RotationOnlySacProblem(adapter));

  // Create a ransac solver for the problem
  opengv::sac::Ransac<RotationOnlySacProblem> ransac;

  ransac.sac_model_ = relposeproblem_ptr;
  ransac.threshold_ = threshold;
  ransac.max_iterations_ = max_iterations;

  // Solve
  ransac.computeModel();
  return arrayFromRotation(ransac.model_coefficients_);
}

} // namespace relative_pose

namespace triangulation
{

bp::object triangulate( bpn::array &b1,
                        bpn::array &b2,
                        bpn::array &t12,
                        bpn::array &R12 )
{
  pyopengv::relative_pose::CentralRelativeAdapter adapter(b1, b2, t12, R12);

  opengv::points_t points;
  for (size_t i = 0; i < adapter.getNumberCorrespondences(); ++i)
  {
    opengv::point_t p = opengv::triangulation::triangulate(adapter, i);
    points.push_back(p);
  }
  return arrayFromPoints(points);
}

bp::object triangulate2( bpn::array &b1,
                         bpn::array &b2,
                         bpn::array &t12,
                         bpn::array &R12 )
{
  pyopengv::relative_pose::CentralRelativeAdapter adapter(b1, b2, t12, R12);

  opengv::points_t points;
  for (size_t i = 0; i < adapter.getNumberCorrespondences(); ++i)
  {
    opengv::point_t p = opengv::triangulation::triangulate2(adapter, i);
    points.push_back(p);
  }
  return arrayFromPoints(points);
}


} // namespace triangulation

} // namespace pyopengv

BOOST_PYTHON_MODULE(pyopengv) {
  using namespace boost::python;

  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  numpy_import_array_wrapper();

  def("absolute_pose_p2p", pyopengv::absolute_pose::p2p);
  def("absolute_pose_p3p_kneip", pyopengv::absolute_pose::p3p_kneip);
  def("absolute_pose_p3p_gao", pyopengv::absolute_pose::p3p_gao);
  def("absolute_pose_gp3p", pyopengv::absolute_pose::gp3p);
  def("absolute_pose_epnp", pyopengv::absolute_pose::epnp);
  def("absolute_pose_gpnp", pyopengv::absolute_pose::gpnp);
  def("absolute_pose_upnp", pyopengv::absolute_pose::upnp);
  def("absolute_pose_optimize_nonlinear", pyopengv::absolute_pose::optimize_nonlinear);
  def("absolute_pose_ransac", pyopengv::absolute_pose::ransac);

  def("relative_pose_twopt", pyopengv::relative_pose::twopt);
  def("relative_pose_twopt_rotation_only", pyopengv::relative_pose::twopt_rotationOnly);
  def("relative_pose_rotation_only", pyopengv::relative_pose::rotationOnly);
  def("relative_pose_fivept_nister", pyopengv::relative_pose::fivept_nister);
  def("relative_pose_fivept_kneip", pyopengv::relative_pose::fivept_kneip);
  def("relative_pose_sevenpt", pyopengv::relative_pose::sevenpt);
  def("relative_pose_eightpt", pyopengv::relative_pose::eightpt);
  def("relative_pose_eigensolver", pyopengv::relative_pose::eigensolver);
  def("relative_pose_sixpt", pyopengv::relative_pose::sixpt);
  def("relative_pose_optimize_nonlinear", pyopengv::relative_pose::optimize_nonlinear);
  def("relative_pose_ransac", pyopengv::relative_pose::ransac);
  def("relative_pose_ransac_rotation_only", pyopengv::relative_pose::ransac_rotationOnly);

  def("triangulation_triangulate", pyopengv::triangulation::triangulate);
  def("triangulation_triangulate2", pyopengv::triangulation::triangulate2);
}
