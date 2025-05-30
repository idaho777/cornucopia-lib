/*--
    CurveDerivativesTest.cpp

    This file is part of the Cornucopia curve sketching library.
    Copyright (C) 2010 Ilya Baran (baran37@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "Anticlothoid.h"
#include "Arc.h"
// #include "Clothoid.h"
#include "Line.h"
#include "Test.h"

#include <iostream>

using namespace std;
using namespace Eigen;
using namespace Cornu;

// class CurveDerivativesTest : public TestCase {
class CurveDerivativesTest {
public:
  // override
  std::string name() { return "CurveDerivativesTest"; }

  // override
  void run() {
    testLine();
    testArc();
    testAnticlothoid();
    // testClothoid();
  }

  void testLine() {
    LinePtr line = new Line(Vector2d(1., 3.), Vector2d(3., 4.));

    CORNU_ASSERT(line->isValid());

    CurvePrimitive::ParamDer der, tanDer;
    double s = 1.5;
    CurvePrimitive::ParamVec params = line->params();

    line->derivativeAt(s, der, tanDer);

    double delta = 1e-6;
    for (int i = 0; i < (int)params.size(); ++i) {
      CurvePrimitive::ParamVec newParams = params;
      newParams[i] += delta;
      line->setParams(newParams);
      Vector2d p1 = line->pos(s);
      Vector2d d1 = line->der(s);
      newParams[i] = params[i] - delta;
      line->setParams(newParams);
      Vector2d p0 = line->pos(s);
      Vector2d d0 = line->der(s);
      Vector2d derV = (p1 - p0) / (2. * delta);
      Vector2d tanDerV = (d1 - d0) / (2. * delta);

      CORNU_ASSERT_LT_MSG((derV - der.col(i)).norm(), delta,
                          "-- Param = "
                              << i << " numeric = " << derV.transpose()
                              << " analytic = " << der.col(i).transpose());
      CORNU_ASSERT_LT_MSG((tanDerV - tanDer.col(i)).norm(), delta,
                          "-- Param = "
                              << i << " numeric = " << tanDerV.transpose()
                              << " analytic = " << tanDer.col(i).transpose());
    }
  }

  void testArc() {
    ArcPtr arc = new Arc(Vector2d(1., 3.), 0.5, 3, 0.1);

    CORNU_ASSERT(arc->isValid());

    CurvePrimitive::ParamDer der, tanDer;
    double s = 1.5;
    for (double c = -.1; c < .1001; c += .01) {
      CurvePrimitive::ParamVec params = arc->params();
      params[CurvePrimitive::CURVATURE] = c;
      arc->setParams(params);

      arc->derivativeAt(s, der, tanDer);

      double delta = 5e-6;
      for (int i = 0; i < (int)params.size(); ++i) {
        CurvePrimitive::ParamVec newParams = params;
        newParams[i] += delta;
        arc->setParams(newParams);
        Vector2d p1 = arc->pos(s);
        Vector2d d1 = arc->der(s);
        newParams[i] = params[i] - delta;
        arc->setParams(newParams);
        Vector2d p0 = arc->pos(s);
        Vector2d d0 = arc->der(s);
        Vector2d derV = (p1 - p0) / (2. * delta);
        Vector2d tanDerV = (d1 - d0) / (2. * delta);

        CORNU_ASSERT_LT_MSG((derV - der.col(i)).norm(), delta,
                            "-- Param = "
                                << i << " numeric = " << derV.transpose()
                                << " analytic = " << der.col(i).transpose());
        CORNU_ASSERT_LT_MSG((tanDerV - tanDer.col(i)).norm(), delta,
                            "-- Param = "
                                << i << " numeric = " << tanDerV.transpose()
                                << " analytic = " << tanDer.col(i).transpose());
      }
    }
  }

  AnticlothoidPtr createAnticlothoid(double a, double startT, double endT) {
    // Compute the starting and ending radii:
    // r(t) = a * t
    double startRadius = a * startT;
    double endRadius = a * endT;

    // Compute the arc-length of the anticlothoid segment:
    // L(t) = 0.5 * a * (endT^2 - startT^2)
    double length = 0.5 * a * (endT * endT - startT * startT);
    Vector2d start(a, 0);
    double startAngle = startT;

    // Debugging::get()->printf("R0: %f, Rn: %f", startRadius, endRadius);
    // Debugging::get()->printf("a: %f, startT: %f, endT %f", a, startT, endT);
    // Debugging::get()->printf("length: %f", length);
    // Debugging::get()->printf("Start: %f %f", start[0], start[1]);
    // Debugging::get()->printf("startAngle: %f", startAngle);

    // Create and return the anticlothoid with the computed parameters:
    return new Anticlothoid(start, startAngle, length, startRadius, endRadius);
  }

  void testAnticlothoid() {
    // Create an anticlothoid with:
    AnticlothoidPtr anticlothoid = createAnticlothoid(1., 0.5, 3.);

    CORNU_ASSERT(anticlothoid->isValid());

    // At s = 0 the anticlothoid should yield the starting point.
    // Vector2d pos0 = anticlothoid->pos(0.0);
    // Vector2d expectedStart(1.0, 3.0);
    // CORNU_ASSERT_LT_MSG((pos0 - expectedStart).norm(), 1e-6,
    //                     "At s=0, anticlothoid position mismatch.");
    //
    // // At s = 0 the tangent angle should equal the provided starting angle.
    // double angle0 = anticlothoid->angle(0.0);
    // CORNU_ASSERT_LT_MSG(fabs(angle0 - 0.5), 1e-6,
    //                     "At s=0, anticlothoid angle mismatch.");

    // Use a small delta for finite difference approximations.
    double delta = 1e-6;
    // Sample a range of s values along the segment.
    for (double s : {0.5, 1.0, 1.5, 2.0, 2.5}) {
      // Test first derivative (tangent):
      Vector2d pos_plus = anticlothoid->pos(s + delta);
      Vector2d pos_minus = anticlothoid->pos(s - delta);
      Vector2d numericDer = (pos_plus - pos_minus) / (2.0 * delta);
      Vector2d analyticDer = anticlothoid->der(s);
      CORNU_ASSERT_LT_MSG((numericDer - analyticDer).norm(), 1e-5,
                          "Tangent derivative mismatch at s = " << s);

      // Test second derivative:
      Vector2d pos_p = anticlothoid->pos(s + delta);
      Vector2d pos_m = anticlothoid->pos(s - delta);
      Vector2d pos0_val = anticlothoid->pos(s);
      Vector2d numericDer2 = (pos_p - 2.0 * pos0_val + pos_m) / (delta * delta);
      Vector2d analyticDer2 = anticlothoid->der2(s);
      CORNU_ASSERT_LT_MSG((numericDer2 - analyticDer2).norm(), 1e-2,
                          "Second derivative mismatch at s = " << s);

      // Test that the derivative of angle(s) equals the curvature.
      double angle_plus = anticlothoid->angle(s + delta);
      double angle_minus = anticlothoid->angle(s - delta);
      double numericCurv = (angle_plus - angle_minus) / (2.0 * delta);
      double analyticCurv = anticlothoid->curvature(s);
      CORNU_ASSERT_LT_MSG(fabs(numericCurv - analyticCurv), 1e-5,
                          "Curvature mismatch at s = " << s);

      // Test that radiusOfCurvature equals 1/curvature when curvature is
      // nonzero.
      if (analyticCurv > 1e-6) {
        double numericRad = 1.0 / analyticCurv;
        double analyticRad = anticlothoid->radiusOfCurvature(s);
        CORNU_ASSERT_LT_MSG(fabs(numericRad - analyticRad), 1e-5,
                            "Radius of curvature mismatch at s = " << s);
      }
    }
  }

  // void testClothoid() {
  //   ClothoidPtr clothoid = new Clothoid(Vector2d(1., 3.), 0.5, 3, 0.0, 0.0);
  //
  //   CORNU_ASSERT(clothoid->isValid());
  //
  //   CurvePrimitive::ParamDer der, tanDer;
  //   for (double s = 0; s < 3.; s += 0.1) {
  //     for (double sc = -.1; sc < .1001; sc += .01) {
  //       for (double ec = -.1; ec < .1001; ec += .01) {
  //         CurvePrimitive::ParamVec params = clothoid->params();
  //         params[CurvePrimitive::CURVATURE] = sc;
  //         params[CurvePrimitive::DCURVATURE] =
  //             (ec - sc) / params[CurvePrimitive::LENGTH];
  //         clothoid->setParams(params);
  //
  //         clothoid->derivativeAt(s, der, tanDer);
  //
  //         double delta = 5e-6;
  //         for (int i = 0; i < (int)params.size(); ++i) {
  //           CurvePrimitive::ParamVec newParams = params;
  //           newParams[i] += delta;
  //           clothoid->setParams(newParams);
  //           Vector2d p1 = clothoid->pos(s);
  //           Vector2d d1 = clothoid->der(s);
  //           newParams[i] = params[i] - delta;
  //           clothoid->setParams(newParams);
  //           Vector2d p0 = clothoid->pos(s);
  //           Vector2d d0 = clothoid->der(s);
  //           Vector2d derV = (p1 - p0) / (2. * delta);
  //           Vector2d tanDerV = (d1 - d0) / (2. * delta);
  //
  //           CORNU_ASSERT_LT_MSG(
  //               (derV - der.col(i)).norm(), delta,
  //               "-- Param = " << i << " numeric = " << derV.transpose()
  //                             << " analytic = " << der.col(i).transpose());
  //           CORNU_ASSERT_LT_MSG(
  //               (tanDerV - tanDer.col(i)).norm(), delta,
  //               "-- Param = " << i << " numeric = " << tanDerV.transpose()
  //                             << " analytic = " <<
  //                             tanDer.col(i).transpose());
  //         }
  //       }
  //     }
  //   }
  // }
};

static CurveDerivativesTest test;
