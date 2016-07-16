//
// Created by kerail on 14.07.16.
//
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>

#include <boost/function.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "geom/SO3.hpp"

typedef boost::multiprecision::cpp_dec_float_100 float_100;
typedef Eigen::Matrix<float_100, 3, 3> Matrix3f100;
typedef SO3<float_100> ReferenceSO3Impl;

typedef std::vector<Eigen::Vector3d> Vector3dList;
typedef std::vector<Matrix3f100> Matrix3f100List;
typedef std::map<std::string, Matrix3f100List> Matrix3f100ListMap;
typedef std::map<std::string, Vector3dList> Vector3dListMap;

using std::cout;
using std::endl;

class AbstractTestSO3 {
 public:
  AbstractTestSO3(const Eigen::Vector3d &v) { }
  virtual Matrix3f100 getMatrix() = 0;
  static std::string name();
};

//class ReferenceSO3: public AbstractTestSO3 {
// public:
//  ReferenceSO3(const Eigen::Vector3d &v) : AbstractTestSO3(v), m_so3(v.cast<float_100>()) { }
//  virtual Matrix3f100 getMatrix() {
//    return m_so3.getQuaternion().toRotationMatrix();
//  }
//  static std::string name() {
//    return "Reference";
//  }
//  ReferenceSO3Impl m_so3;
//};
//
//class ReferenceSO3: public AbstractTestSO3 {
// public:
//  ReferenceSO3(const Eigen::Vector3d &v) : AbstractTestSO3(v) {
//    Eigen::Matrix<float_100, 3, 1> axis = v.cast<float_100>();
//    float_100 angle = axis.norm();
//    if (angle > 0) {
//      axis /= angle;
//    }
//    m_aa = Eigen::AngleAxis<float_100>(angle, axis);
//  }
//  virtual Matrix3f100 getMatrix() {
//    return m_aa.toRotationMatrix();
//  }
//  static std::string name() {
//    return "Reference";
//  }
//  Eigen::AngleAxis<float_100> m_aa;
//};
class ReferenceSO3: public AbstractTestSO3 {
 public:
  ReferenceSO3(const Eigen::Vector3d &v) : AbstractTestSO3(v) {
    Eigen::Matrix<float_100, 3, 1> axis = v.cast<float_100>();
    float_100 angle = axis.norm();
    if (angle > 0) {
      axis /= angle;
    }
    m_aa = Eigen::AngleAxis<float_100>(angle, axis);
  }
  virtual Matrix3f100 getMatrix() {
    return Eigen::Quaternion<float_100>(m_aa).toRotationMatrix();
  }
  static std::string name() {
    return "Reference";
  }
  Eigen::AngleAxis<float_100> m_aa;
};

class TestMySO3Direct: public AbstractTestSO3 {
 public:
  TestMySO3Direct(const Eigen::Vector3d &v) : AbstractTestSO3(v), m_so3(v) { }
  virtual Matrix3f100 getMatrix() {
    return m_so3.getMatrix().cast<float_100>();
  }

  static std::string name() {
    return "My.Direct";
  }
  SO3d m_so3;
};

class TestMySO3Quat: public AbstractTestSO3 {
 public:
  TestMySO3Quat(const Eigen::Vector3d &v) : AbstractTestSO3(v), m_so3(double(v.cast<float_100>().norm()),
                                                                      v.cast<float_100>().norm() > 0
                                                                      ? v.cast<float_100>().normalized().cast<double>()
                                                                      : v) { }
  virtual Matrix3f100 getMatrix() {
    return m_so3.getQuaternion().toRotationMatrix().cast<float_100>();
  }
  static std::string name() {
    return "My.Quater";
  }

  SO3d m_so3;
};

class TestAngleAxisSO3: public AbstractTestSO3 {
 public:
  TestAngleAxisSO3(const Eigen::Vector3d &v)
      : AbstractTestSO3(v), m_aa(v.norm(),
                                 v.norm() > 0 ? v.normalized().cast<double>()
                                              : v) { }
  virtual Matrix3f100 getMatrix() {
    return m_aa.toRotationMatrix().cast<float_100>();
  }

  static std::string name() {
    return "Eigen.AA";
  }
  Eigen::AngleAxisd m_aa;
};

template<size_t TEST_COUNT, class Reference, class... T>
class PrecisionTests;

template<size_t TEST_COUNT, class Reference, class T, class... Ts>
class PrecisionTests<TEST_COUNT, Reference, T, Ts...>: public PrecisionTests<TEST_COUNT, Reference, Ts...> {
  BOOST_STATIC_ASSERT(std::is_base_of<AbstractTestSO3, T>::value);
 public:
  void run() {
    PrecisionTests<TEST_COUNT, Reference, Ts...>::run();
    Vector3dListMap::const_iterator tests_it = getTests().begin();
    const Matrix3f100ListMap &refMap = getReferences();
    for (; tests_it != getTests().end(); ++tests_it) {
      Matrix3f100ListMap::const_iterator refIt = refMap.find(tests_it->first);
      cout << T::name() << "\t" << tests_it->first;
      runTests(tests_it->second, refIt->second);
      cout << endl;
    }
  };
 protected:
  void runTests(const Vector3dList &v, const Matrix3f100List &references) {
    std::vector<float_100> errors(v.size());
#pragma omp parallel for shared(v, errors)
    for (size_t i = 0; i < v.size(); ++i) {
      Eigen::Vector3d coeffs(v[i]);
      T so3(coeffs);
      errors[i] = matrixDiff(references[i], so3.getMatrix());
    }
    float_100 avg, stdDev, stdDevSq, maxErr, minErr;
    calcMetrics(errors, avg, stdDevSq, minErr, maxErr);
    std::cout << "\t" << avg << "\t" << sqrt(stdDevSq)
    << "\t" << minErr << "\t" << maxErr;
  };

  const Vector3dListMap &getTests() {
    return PrecisionTests<TEST_COUNT, Reference, Ts...>::getTests();
  }
  const Matrix3f100ListMap &getReferences() {
    return PrecisionTests<TEST_COUNT, Reference, Ts...>::getReferences();
  }

 private :
  float_100 matrixDiff(const Matrix3f100 &a, const Matrix3f100 &b, const uint8_t type = 2) {
    switch (type) {
      case 0:
        return (a - b).lpNorm<Eigen::Infinity>();
      case 1:
        return (a - b).lpNorm<1>();
      case 2:
      default:
        return (a - b).norm();
    }
  }
  void calcMetrics(const std::vector<float_100> &v,
                   float_100 &avg, float_100 &stdDev,
                   float_100 &minErr, float_100 &maxErr
  ) {
    BOOST_ASSERT(v.size() > 0);
    float_100 sum = 0;
    float_100 divider = 1 / float_100(v.size());
    minErr = v[0];
    maxErr = v[0];
    for (size_t i = 0; i < v.size(); ++i) {
      sum += v[i];
      minErr = std::min<float_100>(v[i], minErr);
      maxErr = std::max<float_100>(v[i], maxErr);
    }
    sum *= divider;
    avg = sum;

    sum = 0;
    for (size_t i = 0; i < v.size(); ++i) {
      float_100 t = (v[i] - avg);
      sum += t * t;
    }
    sum *= divider;
    stdDev = sum;
  }
};

template<size_t TEST_COUNT, class Reference>
class PrecisionTests<TEST_COUNT, Reference> {
  BOOST_STATIC_ASSERT(TEST_COUNT > 0);
 public:
  void run() {
    generateTests();
    calculateReference();
    std::cout << "Method\tTest\tAVG Error\tStdDev Err\tMin Error\tMax Error" << endl;
  };
 protected:

  const Vector3dListMap &getTests() {
    return m_tests;
  }
  const Matrix3f100ListMap &getReferences() {
    return m_referenceResults;
  }

 private:
  void generateTests() {
    cout << "Generation of Tests" << endl;
    m_tests.insert(std::make_pair("USUAL", Vector3dList(TEST_COUNT)));
    m_tests.insert(std::make_pair("ZERO15", Vector3dList(TEST_COUNT)));
    m_tests.insert(std::make_pair("ZERO18", Vector3dList(TEST_COUNT)));
    m_tests.insert(std::make_pair("ZERO22", Vector3dList(TEST_COUNT)));
    Vector3dListMap::iterator usualIt = m_tests.find("USUAL");
    Vector3dListMap::iterator zero15It = m_tests.find("ZERO15");
    Vector3dListMap::iterator zero18It = m_tests.find("ZERO18");
    Vector3dListMap::iterator zero22It = m_tests.find("ZERO22");

#pragma omp parallel for shared(usualIt, zero15It, zero18It, zero22It)
    for (size_t i = 0; i < TEST_COUNT; ++i) {
      Eigen::Vector3d rnd = Eigen::Vector3d::Random();
      rnd.normalize();
      rnd *= M_PI * 2;
      usualIt->second[i] = rnd;
      zero15It->second[i] = rnd * 1e-15;
      zero18It->second[i] = rnd * 1e-18;
      zero22It->second[i] = rnd * 1e-22;
    }
  }

  void calculateReference() {
    cout << "Calculation of precise reference" << endl;
    for (Vector3dListMap::iterator it = m_tests.begin(); it != m_tests.end(); it++) {
      cout << "\tTest: " << it->first << endl;
      m_referenceResults.insert(make_pair(it->first, Matrix3f100List(TEST_COUNT)));
      Matrix3f100ListMap::iterator vIt = m_referenceResults.find(it->first);
#pragma omp parallel for shared(vIt, it)
      for (size_t i = 0; i < TEST_COUNT; ++i) {
        Eigen::Vector3d coeffs(it->second[i]);
        Reference so3Ref(coeffs);
        vIt->second[i] = so3Ref.getMatrix();
      }
    }

  }

  Matrix3f100ListMap m_referenceResults;
  Vector3dListMap m_tests;
};

int main() {
  PrecisionTests<size_t(1000000), ReferenceSO3, TestMySO3Direct, TestMySO3Quat, TestAngleAxisSO3>
      testMachine;
  testMachine.run();
}

