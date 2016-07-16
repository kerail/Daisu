//
// Created by kerail on 10.07.16.
//

#ifndef DAISU_SO3_HPP
#define DAISU_SO3_HPP

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <limits>

template<class T>
class SO3 {
 public:
  typedef Eigen::Matrix<T, 3, 1> Vector3;
  typedef Eigen::Matrix<T, 3, 3> Matrix33;
  typedef Eigen::Quaternion<T> Quaternion;

  SO3(const T &x, const T &y, const T &z) {
    m_coeffs = Eigen::Vector3d(x, y, z);
  }

  SO3(const Vector3 &v) {
    m_coeffs = v;
  }

  SO3(T angle, const Vector3 &v) {
    m_coeffs = v * angle;
  }

  SO3(const Quaternion &q_) {
    const Quaternion &q = q_.w() >= T(0) ? q_ : Quaternion(-q_.coeffs());

    const Vector3 &qv = q.vec();
    T sinha = qv.norm();
    if (sinha > T(0)) {
      T angle = T(2) * atan2(sinha, q.w()); //NOTE: signed
      m_coeffs = qv * (angle / sinha);
    } else {
      // if l is too small, its norm can be equal 0 but norm_inf greater 0
      // probably w is much bigger that vec, use it as length
      m_coeffs = qv * (T(2) / q.w()); ////NOTE: signed
    }
  }

  Vector3 rotateVector(const Vector3 &v) const {
    T a = getAngle();
    Vector3 axis = getAxis();
    T c_a = cos(a);
    T s_a = sin(a);
    return v * c_a + (axis.cross(v)) * s_a + axis * (axis.dot(v)) * (1 - c_a);
  }

  Matrix33 getMatrix() const {
    Matrix33 m;
    const Vector3 &cfs = m_coeffs;

    T l2 = cfs.squaredNorm();
    T l = sqrt(l2);

    T sn_l, cs1_ll, cs;
    if (l == T(0)) {
      // if l is 0 sin(x)/x = 1
      sn_l = T(1);
    } else {
      sn_l = sin(l) / l;
    }

    cs = cos(l);
    static const T c_pi64 = cos(M_PI / T(64));
    if (cs > c_pi64) {//fabs(l) < M_PI/T(64)
      // when l is near nezo, we need to switch to more precise formula
      if (l2 == T(0)) {
        // when l2 is zero, we can precisely calculate limit
        cs1_ll = 1 / T(2);
      } else {
        // 1 - cos(x) = 2 * sin(x/2)^2
        T sn = sin(l / T(2));
        cs1_ll = T(2) * sn * sn / l2;
      }
    } else {
      // here l2 > 0 because abs(l) > pi/64
      cs1_ll = (T(1) - cs) / l2;
    }

    Vector3 sn_ax = sn_l * m_coeffs;
    Vector3 cs1_l_ax = cs1_ll * m_coeffs;

    T tmp;
    tmp = cs1_l_ax.x() * m_coeffs.y();
    m.coeffRef(0, 1) = tmp - sn_ax.z();
    m.coeffRef(1, 0) = tmp + sn_ax.z();

    tmp = cs1_l_ax.x() * m_coeffs.z();
    m.coeffRef(0, 2) = tmp + sn_ax.y();
    m.coeffRef(2, 0) = tmp - sn_ax.y();

    tmp = cs1_l_ax.y() * m_coeffs.z();
    m.coeffRef(1, 2) = tmp - sn_ax.x();
    m.coeffRef(2, 1) = tmp + sn_ax.x();

    m.diagonal() = (cs1_l_ax.cwiseProduct(m_coeffs)).array() + cs;

//    Vector3 sq = cs1_l_ax.cwiseProduct(m_coeffs);
//    m.coeffRef(0, 0) = 1 + (-sq.y() - sq.z());
//    m.coeffRef(1, 1) = 1 + (-sq.x() - sq.z());
//    m.coeffRef(2, 2) = 1 + (-sq.x() - sq.y());

    // Rotation matrix checks
    assert((m * m.transpose() - Matrix33::Identity()).norm() < 1e-10);
    assert(fabs(fabs(m.determinant()) - 1) < 1e-10);
    return m;
  }

  Quaternion getQuaternion() const {
    T a = getAngle();
    Vector3 axis = getAxis();
    T c_a2 = cos(a / 2.0);
    T s_a2 = sin(a / 2.0);
    Vector3 im = axis * s_a2;
    T real = c_a2;

    return Quaternion(real, im.x(), im.y(), im.z());
  }

  Vector3 coeffs() const {
    return m_coeffs;
  }

  T getAngle() const {
    return m_coeffs.norm();
  }

  Vector3 getAxis() const {
    T a = getAngle();
    return a > 0 ? (Vector3) (m_coeffs / a) : m_coeffs;
  }

  SO3<T> inverted() const {
    return SO3<T>(-m_coeffs);
  }

 private:
  Vector3 m_coeffs;

  Matrix33 antiSymmMatrix() const {
    Vector3 axis = getAxis();
    Matrix33 m = Matrix33::Zero();
    m(0, 1) = -axis(2);
    m(1, 0) = axis(2);

    m(0, 2) = axis(1);
    m(2, 0) = -axis(1);

    m(2, 1) = axis(0);
    m(1, 2) = -axis(0);
    return m;
  }

};

typedef SO3<double> SO3d;
typedef SO3<float> SO3f;

#endif //DAISU_SO3_HPP
