//
// Created by kerail on 11.07.16.
//

#include "geom/SO3.hpp"

#include <boost/function.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include <iostream>

typedef boost::multiprecision::cpp_dec_float_100 float_100;
typedef Eigen::Matrix<float_100, 3, 1> Vector3f100;
typedef Eigen::Matrix<float_100, 3, 3> Matrix3f100;
typedef Eigen::Quaternion<float_100> Quaternionf100;
typedef SO3<float_100> SO3f100;

inline void so3_test_matrix(const SO3d &s, const Eigen::Vector3d &v) {
  Eigen::Vector3d v1 = s.getMatrix() * v;
}

inline void so3_test_quaternion(const SO3d &s, const Eigen::Vector3d &v) {
  Eigen::Vector3d v1 = s.getQuaternion().toRotationMatrix() * v;
}

inline void so3_test_rotation(const SO3d &s, const Eigen::Vector3d &v) {
  Eigen::Vector3d v1 = s.rotateVector(v);
}

double test(size_t testCount, const boost::function<void(const SO3d &, const Eigen::Vector3d &)> &f) {
  clock_t startClock = clock();
  for (size_t i = 0; i < testCount; ++i) {
    SO3d so3((Eigen::Vector3d) Eigen::Vector3d::Random());
    Eigen::Vector3d testVec = Eigen::Vector3d::Random();

    so3_test_matrix(so3, testVec);
  }

  clock_t endClock = clock();
  return (endClock - startClock) / double(CLOCKS_PER_SEC);
}

int main() {
  const size_t testCount = 1e8;

  std::cout << "SO3 matrix time: " << test(testCount, &so3_test_matrix) << std::endl;
  std::cout << "SO3 quaternion matrix time: " << test(testCount, &so3_test_quaternion) << std::endl;
  std::cout << "SO3 direct rotation time: " << test(testCount, &so3_test_rotation) << std::endl;
}