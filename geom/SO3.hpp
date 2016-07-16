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
    fromQuaternion(q_);
  }

  SO3(const Matrix33 &m) {
    assert(isRotationMatrix(m));
    //TODO: implement
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
    static const T c_pi4096 = cos(M_PI / T(4096));
    if (cs > c_pi4096) {//fabs(l) < M_PI/T(4096)
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
      // here l2 > 0 because abs(l) > pi/4096
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
    assert(isRotationMatrix(m));
    return m;
  }

  Quaternion getQuaternion() const {
    const Vector3 cf2 = m_coeffs / T(2);
    T a = cf2.norm();
    if (a > T(0)) {
      T sn = sin(a) / a;
      return Quaternion(cos(a), cf2.x() * sn, cf2.y() * sn, cf2.z() * sn);
    } else {
      return Quaternion(T(1), cf2.x(), cf2.y(), cf2.z());
    }
  }

  Vector3 rotateVector(const Vector3 &v) const {
    T l2 = m_coeffs.squaredNorm();
    T l = sqrt(l2);

    T sa_l, ca, ca1_ll;
    if (l == T(0)) {
      sa_l = 1;
    } else {
      sa_l = sin(l) / l;
    }

    ca = cos(l);

    static const T c_pi4096 = cos(M_PI / T(4096));
    if (ca > c_pi4096) {//fabs(l) < M_PI/T(4096)
      // when l is near nezo, we need to switch to more precise formula
      if (l2 == T(0)) {
        // when l2 is zero, we can precisely calculate limit
        ca1_ll = 1 / T(2);
      } else {
        // 1 - cos(x) = 2 * sin(x/2)^2
        T sn = sin(l / T(2));
        ca1_ll = T(2) * sn * sn / l2;
      }
    } else {
      // here l2 > 0 because abs(l) > pi/4096
      ca1_ll = (T(1) - ca) / l2;
    }

    return v * ca + (m_coeffs * sa_l).cross(v) + m_coeffs * ((m_coeffs.dot(v)) * ca1_ll);
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

  SO3<T> operator *(const SO3<T> &r) {
    SO3<T> l(*this);
    return l *= r;
  }

  SO3<T> &operator *=(const SO3<T> &l) {
    fromQuaternion(getQuaternion() * l.getQuaternion());
    return *this;
  }

  SO3<T> operator ~() {
    return inverted();
  }

 private:
  Vector3 m_coeffs;

  void fromQuaternion(const Quaternion &q_) {
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
  bool isRotationMatrix(const Matrix33 &m) const {
    return ((m * m.transpose() - Matrix33::Identity()).norm() < 1e-10)
           && (fabs(fabs(m.determinant()) - 1) < 1e-10);
  }
};

typedef SO3<double> SO3d;
typedef SO3<float> SO3f;

#endif //DAISU_SO3_HPP
