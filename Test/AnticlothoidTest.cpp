/*--
    AnticlothoidTest.cpp

    Copyright (C) 2025 Joonho Kim (idaho777@gmail.com)

    Previous documentation:
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

#include "Debugging.h"
#include "Test.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "Anticlothoid.h"
#include "PrimitiveFitUtils.h"

using namespace std;
using namespace Eigen;
using namespace Cornu;

class AnticlothoidTest : public TestCase {
  enum Param { X = 0, Y, ANGLE, LENGTH, RADIUS, ARADIUS };

  struct ACP {
    double a;
    double startT;
    double endT;
    double noise;
    Vector2d T;
    double theta;
    bool flip;

    ACP(double a, double startT, double endT, double noise, const Vector2d &T,
        double theta, bool flip)
        : a(a), startT(startT), endT(endT), noise(noise), T(T), theta(theta),
          flip(flip) {}
  };

public:
  // override
  std::string name() { return "Anticlothoid Super Test"; }

  // override
  void run() {
    vector<ACP> params = createAnticlothoids();

    for (const ACP &param : params) {
      // Debugging::get()->printf(
      //     "ANTICLOTHOIDTEST: a:%6.2f, startT:%6.2f, endT:%6.2f, noise:%6.2f,
      //     " "Translate:(%6.2f, %6.2f), rotate:%6.2f, flip: %d", param.a,
      //     param.startT, param.endT, param.noise, param.T[0], param.T[1],
      //     param.theta, param.flip);

      // Test the anticlothoid against the sampled points:
      // testAnticlothoidStart(param);
      // testAnticlothoidDerivatives(param);
      // testAnticlothoidAngleCurvature(param);
      // testAnticlothoidProject(param);
      // testAnticlothoidTrim(param);
      // testAnticlothoidFlip(param);
      testAnticlothoidFitter(param);
    }
  }

  vector<ACP> createAnticlothoids() {
    vector<ACP> params;
    // a, startT, endT, noise, Translate, rotate, flip

    // Canonical Anticlothoid
    vector<Vector2d> Ts;
    Ts.push_back(Vector2d(0.0, 0.0));
    Ts.push_back(Vector2d(2.0, 3.0));
    Ts.push_back(Vector2d(4.0, -5.0));
    Ts.push_back(Vector2d(-10.0, -2.0));
    Ts.push_back(Vector2d(-8.0, 3.0));

    vector<double> thetas;
    thetas.push_back(0.0);
    thetas.push_back(M_PI);
    thetas.push_back(M_PI / 4.0);
    thetas.push_back(2 * M_PI);
    thetas.push_back(M_PI / 3);

    vector<double> as;
    as.push_back(1.0);
    as.push_back(0.4);
    as.push_back(2);
    as.push_back(4);
    as.push_back(16);

    for (double a : as) {
      for (Vector2d T : Ts) {
        for (double theta : thetas) {
          for (int i = 1; i >= -1; i -= 2) {
            double neg = i * 1.0;
            params.emplace_back(a, 0.0, 3.0, 0.0, T, neg * theta, false);
            params.emplace_back(a, 2.0, 5.0, 0.0, T, neg * theta, false);
            params.emplace_back(a, 1, 10, 0.0, T, neg * theta, false);
            params.emplace_back(a, 3, 8, 0.0, T, neg * theta, false);
          }
        }
      }
    }

    return params;
  }

  AnticlothoidPtr createAnticlothoid(ACP params) {
    double a = params.a, startT = params.startT, endT = params.endT;
    // Compute the starting and ending radii:
    // r(t) = a * t
    double startRadius = a * startT;
    double endRadius = a * endT;

    // Compute the arc-length of the anticlothoid segment:
    // L(t) = 0.5 * a * (endT^2 - startT^2)
    double length = 0.5 * a * (endT * endT - startT * startT);
    Vector2d start;
    start[0] = a * (std::cos(startT) + startT * std::sin(startT));
    start[1] = a * (std::sin(startT) - startT * std::cos(startT));

    double startAngle = startT;

    Vector2d T = params.T;
    double theta = params.theta;

    Matrix2d rot;
    rot << cos(theta), -sin(theta), sin(theta), cos(theta);
    start = rot * start + T;

    startAngle += theta;

    // Debugging::get()->printf("R0: %f, Rn: %f", startRadius, endRadius);
    // Debugging::get()->printf("a: %f, startT: %f, endT %f", a, startT, endT);
    // Debugging::get()->printf("length: %f", length);
    // Debugging::get()->printf("Start: %f %f", start[0], start[1]);
    // Debugging::get()->printf("startAngle: %f", startAngle);

    // Create and return the anticlothoid with the computed parameters:
    return new Anticlothoid(start, startAngle, length, startRadius, a);
  }

  // vector<Vector2d> sample_points(double a, double startT, double endT,
  //                                double noise) {
  //   AnticlothoidPtr ac = createAnticlothoid(a, startT, endT);
  //   CORNU_ASSERT(ac->isValid());
  //
  //   vector<Vector2d> points;
  //   for (double s = 0; s < ac->length(); s += 0.1) {
  //     Vector2d p = ac->pos(s);
  //     p[0] += noise * ((double)rand() / (RAND_MAX));
  //     p[1] += noise * ((double)rand() / (RAND_MAX));
  //     points.push_back(p);
  //   }
  //   return points;
  // }

  void getActualPoint(ACP params, double s, Vector2d &p, double &angle) {
    double a = params.a;
    double theta = params.theta;

    double t0 = params.startT;
    double s0 = 0.5 * a * (t0 * t0); // arc-length at t0

    double trueS = s0 + s;

    bool neg = trueS < 0;
    if (neg) {
      trueS *= -1;
      // theta += M_PI;
    }

    double t = sqrt(2 * trueS / a);
    if (neg) {
      t = -t;
    }
    double x = a * (cos(t) + t * sin(t));
    double y = a * (sin(t) - t * cos(t));
    Matrix2d rot;
    rot << cos(theta), -sin(theta), sin(theta), cos(theta);
    Vector2d translate = params.T;

    p = rot * Vector2d(x, y) + translate;
    angle = AngleUtils::toRange(theta + t);
  }

  //------------------------------------------------------------------------------
  // Test 1: Starting Conditions
  //------------------------------------------------------------------------------
  void testAnticlothoidStart(ACP params) {
    AnticlothoidPtr ac = createAnticlothoid(params);
    CORNU_ASSERT(ac->isValid());

    // printAnticlothoid(ac);

    double startS = 0.5 * params.a * (params.startT * params.startT);
    double endS = 0.5 * params.a * (params.endT * params.endT);

    Vector2d pos1;
    double ang1;

    Vector2d actualPos = ac->pos(startS);
    double actualAng = ac->angle(startS);
    getActualPoint(params, startS, pos1, ang1);
    CORNU_ASSERT_LT_MSG((actualPos - pos1).norm(), 1e-6,
                        "Anticlothoid start position mismatch: pos():"
                            << actualPos[0] << " " << actualPos[1]
                            << "get():" << pos1[0] << " " << pos1[1]);
    CORNU_ASSERT_LT_MSG(fabs(actualAng - ang1), 1e-6,
                        "Anticlothoid start angle mismatch: ang():"
                            << actualAng << "get():" << ang1);

    actualPos = ac->pos(endS);
    actualAng = ac->angle(endS);
    getActualPoint(params, endS, pos1, ang1);
    // Debugging::get()->printf("Anticlothoid Actual: %f %f", actualPos[0],
    //                          actualPos[1]);
    // Debugging::get()->printf("Anticlothoid s: %f %f", pos1[0], pos1[1]);
    CORNU_ASSERT_LT_MSG((actualPos - pos1).norm(), 1e-6,
                        "Anticlothoid end position mismatch");
    CORNU_ASSERT_LT_MSG(fabs(actualAng - ang1), 1e-6,
                        "Anticlothoid end angle mismatch");

    actualPos = ac->pos(-endS);
    actualAng = ac->angle(-endS);
    getActualPoint(params, -endS, pos1, ang1);
    // Debugging::get()->printf("Anticlothoid negative Actual: %f %f",
    //                          actualPos[0], actualPos[1]);
    // Debugging::get()->printf("Anticlothoid negative s: %f %f", pos1[0],
    //                          pos1[1]);
    CORNU_ASSERT_LT_MSG((actualPos - pos1).norm(), 1e-6,
                        "Anticlothoid negative position mismatch");
    CORNU_ASSERT_LT_MSG(fabs(actualAng - ang1), 1e-6,
                        "Anticlothoid negative angle mismatch");
  }

  //------------------------------------------------------------------------------
  // Test 2: Derivative Verification
  //------------------------------------------------------------------------------
  void testAnticlothoidDerivatives(ACP params) {
    AnticlothoidPtr ac = createAnticlothoid(params);
    CORNU_ASSERT(ac->isValid());

    double delta = 5e-6;
    // Test at several sample s values.
    for (double s : {0.5, 1.0, 1.5, 2.0, 2.5}) {
      // First derivative test (tangent):
      Vector2d p_plus = ac->pos(s + delta);
      Vector2d p_minus = ac->pos(s - delta);
      Vector2d numericDer = (p_plus - p_minus) / (2.0 * delta);
      Vector2d analyticDer = ac->der(s);
      CORNU_ASSERT_LT_MSG((numericDer - analyticDer).norm(), 1e-5,
                          "First derivative mismatch at s=" << s);

      // Second derivative test:
      Vector2d p1 = ac->pos(s + delta);
      Vector2d p0 = ac->pos(s);
      Vector2d p_1 = ac->pos(s - delta);
      Vector2d numericDer2 = (p1 - 2.0 * p0 + p_1) / (delta * delta);
      Vector2d analyticDer2 = ac->der2(s);
      CORNU_ASSERT_LT_MSG((numericDer2 - analyticDer2).norm(), 1e-3,
                          "Second derivative mismatch at s=" << s);
    }
  }

  //------------------------------------------------------------------------------
  // Test 3: Angle, Curvature, and Radius of Curvature
  //------------------------------------------------------------------------------
  void testAnticlothoidAngleCurvature(ACP params) {
    AnticlothoidPtr ac = createAnticlothoid(params);
    CORNU_ASSERT(ac->isValid());

    double delta = 1e-6;
    // Test at several s values.
    for (double s : {0.5, 1.0, 1.5, 2.0, 2.5}) {
      double a0 = ac->angle(s);
      // Finite difference of angle:
      double a_plus = ac->angle(s + delta);
      double a_minus = ac->angle(s - delta);
      double numericCurv = (a_plus - a_minus) / (2.0 * delta);
      double analyticCurv = ac->curvature(s);
      CORNU_ASSERT_LT_MSG(fabs(numericCurv - analyticCurv), 1e-5,
                          "Curvature mismatch at s=" << s);

      // Test radius of curvature (reciprocal of curvature when curvature != 0)
      if (analyticCurv > 1e-6) {
        double numericRoC = 1.0 / analyticCurv;
        double analyticRoC = ac->radiusOfCurvature(s);
        CORNU_ASSERT_LT_MSG(fabs(numericRoC - analyticRoC), 1e-5,
                            "Radius of curvature mismatch at s=" << s);
      }
    }
  }

  //------------------------------------------------------------------------------
  // Test 4: Projection Function
  //------------------------------------------------------------------------------
  void testAnticlothoidProject(ACP params) {
    AnticlothoidPtr ac = createAnticlothoid(params);
    CORNU_ASSERT(ac->isValid());

    double s_sample = ac->length() / 2;

    // Choose a sample s value and compute its position.
    Vector2d pos = ac->pos(s_sample);
    // Perturb the point a little.
    Vector2d queryPoint = pos + Vector2d(0.05, -0.03);
    double projected_s = ac->project(queryPoint);
    Vector2d projPos = ac->pos(projected_s);
    double dist = (projPos - queryPoint).norm();
    // Check that the distance is small.
    CORNU_ASSERT_LT_MSG(dist, 0.1,
                        "Projection did not return close enough s value");
  }

  //------------------------------------------------------------------------------
  // Test 5: Trim Functionality
  //------------------------------------------------------------------------------
  void testAnticlothoidTrim(ACP params) {
    AnticlothoidPtr ac = createAnticlothoid(params);
    CORNU_ASSERT(ac->isValid());
    double length = ac->length();

    // double new_from = (((double)rand() / (RAND_MAX)) / 2.) * length;
    double new_from = 1;
    // double new_to = (((double)rand() / (RAND_MAX)) / 2. + 0.5) * length;
    double new_to = 3;

    double expectedRadiusFrom = ac->radiusOfCurvature(new_from);
    double expectedRadiusTo = ac->radiusOfCurvature(new_to);

    double expectedAngleFrom = ac->angle(new_from);
    double expectedAngleTo = ac->angle(new_to);

    Vector2d expectedFrom = ac->pos(new_from);
    Vector2d expectedTo = ac->pos(new_to);

    // auto p = ac->params();
    // double _a = ac->getA();
    // Debugging::get()->printf("Params: %f %f %f %f %f %f %f", p[0], p[1],
    // p[2],
    //                          p[3], p[4], p[5], _a);

    ac->trim(new_from, new_to);

    // At s=0 of the trimmed curve, the position should equal the oldpos(s=1.0)
    Vector2d newStart = ac->pos(0.0);
    CORNU_ASSERT_LT_MSG((newStart - expectedFrom).norm(), 1e-6,
                        "Trim: starting point did not update correctly");

    double newLength = ac->length();
    CORNU_ASSERT_LT_MSG(fabs(newLength - (new_to - new_from)), 1e-6,
                        "Trim: LENGTH parameter mismatch");

    Vector2d newEnd = ac->pos(ac->length());
    CORNU_ASSERT_LT_MSG((newEnd - expectedTo).norm(), 1e-6,
                        "Trim: ending point did not update correctly");

    double newAngleStart = ac->angle(0.0);
    CORNU_ASSERT_LT_MSG(fabs(newAngleStart - expectedAngleFrom), 1e-6,
                        "Trim: starting angle mismatch");
    double newAngleEnd = ac->angle(ac->length());
    CORNU_ASSERT_LT_MSG(fabs(newAngleEnd - expectedAngleTo), 1e-6,
                        "Trim: ending angle mismatch");

    double newRadiusStart = ac->radiusOfCurvature(0.0);
    CORNU_ASSERT_LT_MSG(fabs(newRadiusStart - expectedRadiusFrom), 1e-6,
                        "Trim: RADIUS start parameter mismatch");
    double newRadiusEnd = ac->radiusOfCurvature(ac->length());
    CORNU_ASSERT_LT_MSG(fabs(newRadiusEnd - expectedRadiusTo), 1e-6,
                        "Trim: RADIUS end parameter mismatch");
  }

  //------------------------------------------------------------------------------
  // Test 6: Flip Functionality
  //------------------------------------------------------------------------------
  void testAnticlothoidFlip(ACP params) {
#define FLIPDEBUG 0
    AnticlothoidPtr ac = createAnticlothoid(params);
    CORNU_ASSERT(ac->isValid());

    double l = ac->length();

    Vector2d oldStart = ac->pos(0.0);
    Vector2d oldEnd = ac->pos(l);
    double oldStartAngle = ac->angle(0.0);
    double oldEndAngle = ac->angle(l);
    double oldParamR0 = ac->radiusOfCurvature(0.0);
    double oldLength = ac->length();
    double oldStartR = ac->angle(0);
    double oldEndR = ac->radiusOfCurvature(l);

    // Flip the anticlothoid
#if FLIPDEBUG
    Debugging::get()->printf("Before flip:");
    printAnticlothoid(ac);
    Debugging::get()->printf("Flipping...");
#endif
    ac->flip();
#if FLIPDEBUG
    Debugging::get()->printf("After flip:");
    printAnticlothoid(ac);
#endif

    // After flipping, the new starting point should equal the old end point
    // and the new ending point should equal the old start point.
    // A  -> B
    // B' <- A'
    Vector2d newStart = ac->pos(0.0);
    CORNU_ASSERT_LT_MSG((newStart - oldEnd).norm(), 1e-6,
                        "Flip: new starting point mismatch");
    Vector2d newEnd = ac->pos(l);
    CORNU_ASSERT_LT_MSG((newEnd - oldStart).norm(), 1e-6,
                        "Flip: new end point mismatch");

    // newStartAngle should be in the opposite direction of the old end angle
    // and the newEndAngle should be in the opposite direciton of the old start
    // angle.
    double newStartAngle = ac->params()[ANGLE];
    double expectedAngle = AngleUtils::toRange(oldEndAngle + M_PI);

    CORNU_ASSERT_LT_MSG(fabs(newStartAngle - expectedAngle), 1e-6,
                        "Flip: ANGLE mismatch" << newStartAngle << " "
                                               << oldEndAngle);
    // Debugging::get()->printf("New start angle: %f, oldEndAngle: %f",
    //                          newStartAngle, oldEndAngle);
    newStartAngle = ac->angle(0);
    expectedAngle = AngleUtils::toRange(oldEndAngle + M_PI);
    CORNU_ASSERT_LT_MSG(fabs(newStartAngle - expectedAngle), 1e-6,
                        "Flip: angle(l=0) mismatch" << newStartAngle << " "
                                                    << oldEndAngle);

    // Debugging::get()->printf("New end angle: %f, oldStartAngle: %f",
    //                          newEndAngle, oldStartAngle);
    double newEndAngle = ac->angle(l);
    expectedAngle = AngleUtils::toRange(oldStartAngle + M_PI);
    CORNU_ASSERT_LT_MSG(fabs(newEndAngle - expectedAngle), 1e-6,
                        "Flip: ending angle mismatch" << newEndAngle << " "
                                                      << oldStartAngle);

    double length = ac->params()[LENGTH];
    CORNU_ASSERT_LT_MSG(fabs(length - oldLength), 1e-6,
                        "Flip: LENGTH parameter mismatch");

    double r0 = ac->params()[RADIUS];
    CORNU_ASSERT_LT_MSG(
        r0 + oldEndR, 1e-6,
        "Flip: starting radius mismatch r0:" << r0 << " oldEndR:" << oldEndR);

    // Unflip the anticlothoid.
#if FLIPDEBUG
    Debugging::get()->printf("Flipping...");
#endif
    ac->flip();
#if FLIPDEBUG
    Debugging::get()->printf("After Un-flip:");
    printAnticlothoid(ac);
#endif

    newStart = ac->pos(0.0);
    CORNU_ASSERT_LT_MSG((newStart - oldStart).norm(), 1e-6,
                        "UnFlip: new starting point mismatch");
    newEnd = ac->pos(l);
    CORNU_ASSERT_LT_MSG((newEnd - oldEnd).norm(), 1e-6,
                        "UnFlip: new end point mismatch");

    newStartAngle = ac->params()[ANGLE];
    CORNU_ASSERT_LT_MSG(fabs(newStartAngle - oldStartAngle), 1e-6,
                        "UnFlip: starting angle mismatch");

    length = ac->params()[LENGTH];
    CORNU_ASSERT_LT_MSG(fabs(length - oldLength), 1e-6,
                        "UnFlip: LENGTH parameter mismatch");

    r0 = ac->params()[RADIUS];
    CORNU_ASSERT_LT_MSG(fabs(r0 - oldParamR0), 1e-6,
                        "UnFlip: starting radius mismatch");
  }

  //------------------------------------------------------------------------------
  // Test 7: AnticlothoidFitter
  //------------------------------------------------------------------------------
  void testAnticlothoidFitter(ACP params) {
    // Create a set of sample points that roughly lie on an anticlothoid curve.
    double a = params.a;
    double ts[] = {0, 1, 2, 3};
    vector<Vector2d> pts;
    for (double t : ts) {
      double x = a * (cos(t) + t * sin(t));
      double y = a * (sin(t) - t * cos(t));
      pts.emplace_back(x, y);
    }

    // Create an AnticlothoidFitter instance and add points.
    AnticlothoidFitter fitter;
    for (const auto &p : pts) {
      fitter.addPoint(p);
    }

    double candidateError;
    bool pass = fitter.verifyCandidate(candidateError);
    // The error should be non-negative. (Depending on data, error may not be
    // very small.)
    // Debugging::get()->printf("Fitter: candidate error: %f", candidateError);
    // CORNU_ASSERT_MSG(candidateError >= 0,
    //                  "Error is too large: " << candidateError);
    //
    // // Get the fitted anticlothoid.
    // AnticlothoidPtr ac = fitter.getCurve();
    // CORNU_ASSERT(ac->isValid());
    //
    // // Optionally, we can check that the starting point and length roughly
    // match
    // // our data.
    // Vector2d fitStart = ac->pos(0.0);
    // CORNU_ASSERT_LT_MSG((fitStart - pts.front()).norm(), 0.5,
    //                     "Fitter: Starting point is off.");
    // double fitLength = ac->params()[LENGTH];
    // CORNU_ASSERT(fitLength > 0);
  }

  // Test fitting for line
  // Test fitting for Arc
  // Test fitting for Anticlothoid

  void printAnticlothoid(AnticlothoidPtr ac) {
    // Print the test results.
    auto p = ac->params();
    Vector2d shift;
    Matrix2d mat;
    (void)ac->getTransform(mat, shift);
    Debugging::get()->printf(
        "Param:\t v:(% 6.2f, % 6.2f), Theta:% 6.2f, currL:% "
        "6.2f, RoC0:% 6.2f A:% 6.2f, Tr:(% 6.2f, % 6.2f), "
        "Rot:(% 6.2f, % 6.2f, % 6.2f, % 6.2f)",
        p[0], p[1], p[2], p[3], p[4], p[5], shift[0], shift[1], mat(0, 0),
        mat(0, 1), mat(1, 0), mat(1, 1));

    printAnticlothoid(ac, 0);
    printAnticlothoid(ac, ac->length());
  }

  void printAnticlothoid(AnticlothoidPtr ac, const double s) {
    Vector2d p = ac->pos(s);
    double r = ac->radiusOfCurvature(s);
    double theta = ac->angle(s);
    double trueS, t;
    ac->getSandT(s, trueS, t);
    Debugging::get()->printf(
        "    s: %4.1f,\t v:(% 6.2f, % 6.2f), theta:% 6.2f, trueL:% "
        "6.2f, RoCL:% 6.2f t:% 6.2f",
        s, p[0], p[1], theta, trueS, r, t);
  }
};

static AnticlothoidTest test;
