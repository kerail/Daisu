//
// Created by kerail on 10.07.16.
//
#include "geom/SO3.hpp"

#define BOOST_TEST_MODULE boost_test_so3
#define BOOST_TEST_LOG_LEVEL all
#include <boost/test/unit_test.hpp>

using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::AngleAxisd;

BOOST_AUTO_TEST_CASE(test_so3) {
  const double tol = 1e-5;
  boost::unit_test::unit_test_log.set_threshold_level(boost::unit_test::log_messages);

  for (double i = -1; i <= 1; i += 1) {
    for (double j = -1; j <= 1; j += 1) {
      for (double k = -1; k <= 1; k += 1) {
        for (double a = -2; a <= 2; a += 0.5) {
          Vector3d axisAngle(i, j, k);
          if (fabs(i) + fabs(j) + fabs(k) > 0) {
            axisAngle.normalize();
          }
          double angle = M_PI * a;

          BOOST_TEST_MESSAGE("Axis:\n" << axisAngle << "\n");
          BOOST_TEST_MESSAGE("Angle:\n" << angle << "\n");

          BOOST_TEST_MESSAGE(
              "Rotation SO3:\n" << axisAngle * angle << " \nNorm: " << (axisAngle * angle).norm() << "\n");

          SO3d so3(axisAngle * angle);

          Quaterniond q(AngleAxisd(angle, axisAngle));
          q.normalize();

          SO3d so3FromQ(q);

          Quaterniond rotationQ = so3.getQuaternion();

          Quaterniond rotationQinv = so3.inverted().getQuaternion();
          Quaterniond qInv = q.inverse();

          BOOST_TEST_MESSAGE("Rotation from so3:\n" << rotationQ.coeffs() << "\n");
          BOOST_TEST_MESSAGE("Rotation from q:\n" << q.coeffs() << "\n");

          bool t1 = (rotationQ.coeffs() - q.coeffs()).norm() < tol;
          bool t2 = (rotationQ.coeffs() + q.coeffs()).norm() < tol;

          BOOST_REQUIRE(t1 || t2);

          bool t3 = (rotationQinv.coeffs() - qInv.coeffs()).norm() < tol;
          bool t4 = (rotationQinv.coeffs() + qInv.coeffs()).norm() < tol;

          BOOST_REQUIRE(t3 || t4);

          BOOST_TEST_MESSAGE("SO3:\n" << so3.coeffs() << "\n");
          BOOST_TEST_MESSAGE("SO3 from q:\n" << so3FromQ.coeffs() << "\n");

//          BOOST_REQUIRE((so3FromQ.coeffs() - so3.coeffs()).norm() < tol);
          bool t5 = (so3FromQ.getQuaternion().coeffs() - so3.getQuaternion().coeffs()).norm() < tol;
          bool t6 = (so3FromQ.getQuaternion().coeffs() + so3.getQuaternion().coeffs()).norm() < tol;
          BOOST_REQUIRE(t5 || t6);

          for (double i2 = -1; i2 <= 1; i2 += 1) {
            for (double j2 = -1; j2 <= 1; j2 += 1) {
              for (double k2 = -1; k2 <= 1; k2 += 1) {
                Vector3d v(i2, j2, k2);

                Vector3d checkVector = q.toRotationMatrix() * v;

                Vector3d rotatedVecByFunc = so3.rotateVector(v);
                Vector3d rotatedVecByMatrix = so3.getMatrix() * v;

                BOOST_TEST_MESSAGE("Check vector:\n " << checkVector << "\n\n");
                BOOST_TEST_MESSAGE("Check rotatedVecByFunc:\n " << rotatedVecByFunc << "\n");
                BOOST_TEST_MESSAGE("Check rotationMatrix so3:\n " << so3.getMatrix() << "\n");
                BOOST_TEST_MESSAGE("Check rotationMatrix q:\n " << q.toRotationMatrix() << "\n");
                BOOST_TEST_MESSAGE("Check rotatedVecByMatrix:\n " << rotatedVecByMatrix << "\n");

                BOOST_REQUIRE((rotatedVecByFunc - checkVector).norm() < tol);
                BOOST_REQUIRE((rotatedVecByMatrix - checkVector).norm() < tol);

              }
            }
          }
        }
      }
    }
  }


//  cout << "Original vector: " << v << endl << endl;
//  cout << "Rotate vector by func: " << endl << so3.rotateVector(v) << endl << endl;
//  cout << "Rotation matrix: " << endl << so3.getMatrix() << endl << endl;
//  cout << "Rotate vector by matrix: " << endl << (so3.getMatrix()*v) << endl << endl;
//  cout << "Rotate vector by quaternion: " << endl << (so3.getQuaternion().toRotationMatrix()*v) << endl << endl;
//  cout << "Rotate vector by test quaternion: " << endl << (q.toRotationMatrix()*v) << endl << endl;
}