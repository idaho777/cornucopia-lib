/*--
    PrimitiveFitUtils.cpp

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

#include "PrimitiveFitUtils.h"

#include "AngleUtils.h"
#include "Anticlothoid.h"
#include "Arc.h"
#include "Clothoid.h"
#include "Fresnel.h"
#include "Line.h"
#include "VectorC.h"

#include "Debugging.h"

#include <Eigen/Eigenvalues>
#include <cassert>
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

    void
    LineFitter::addPointW(const Vector2d &pt, double weight) {
  ++_numPts;
  _totWeight += weight;

  _lastPoint = pt;
  if (_numPts == 1)
    _firstPoint = pt;

  _sum += weight * pt;
  _squaredSum[0] += weight * SQR(pt[0]);
  _squaredSum[1] += weight * SQR(pt[1]);
  _crossSum += weight * pt[0] * pt[1];
}

LinePtr LineFitter::getCurve() const {
  if (_numPts < 2)
    return LinePtr();

  Matrix2d cov = Matrix2d::Zero();
  double factor = 1. / double(_totWeight);
  cov(0, 0) = _squaredSum[0] - SQR(_sum[0]) * factor;
  cov(1, 1) = _squaredSum[1] - SQR(_sum[1]) * factor;
  cov(0, 1) = cov(1, 0) = _crossSum - _sum[0] * _sum[1] * factor;
  cov *= factor;

  SelfAdjointEigenSolver<Matrix2d> eigenSolver(cov);

  Vector2d eigVs = eigenSolver.eigenvalues();
  Vector2d dir = eigenSolver.eigenvectors()
                     .col(1)
                     .normalized(); // 1 is the index of the larger eigenvalue
  Vector2d pt = _sum * factor;
  Vector2d pt0 = pt + ((_firstPoint - pt).dot(dir)) * dir;
  Vector2d pt1 = pt + ((_lastPoint - pt).dot(dir)) * dir;

  return new Line(pt0, pt1);
}

void ArcFitter::addPointW(const Vector2d &pt, double weight) {
  _pts.push_back(pt);

  Vector3d pt3(pt[0] - _pts[0][0], pt[1] - _pts[0][1],
               (pt - _pts[0]).squaredNorm());

  _totWeight += weight;
  _sum += weight * pt3;
  _squaredSum += weight * pt3 * pt3.transpose();
}

ArcPtr ArcFitter::getCurve() const {
  if ((int)_pts.size() < 2)
    return ArcPtr();

  double factor = 1. / _totWeight;
  Vector3d pt = _sum * factor;
  Matrix3d cov = factor * _squaredSum - pt * pt.transpose();

  SelfAdjointEigenSolver<Matrix3d> eigenSolver(cov);
  Vector3d eigVs = eigenSolver.eigenvalues();

  Vector3d dir = eigenSolver.eigenvectors().col(
      0); // 0 is the index of the smallest eigenvalue
  dir /= (1e-16 + dir[2]);

  double dot = dir.dot(pt);
  // circle equation is:
  // dir[0] * x + dir[1] * y + (x^2+y^2) = dot
  Vector2d center = -0.5 * Vector2d(dir[0], dir[1]);
  double radius = sqrt(1e-16 + dot + center.squaredNorm());
  center += _pts[0];

  // TODO: convert code to use AngleUtils
  // Now get the arc
  Vector2d c[3] = {_pts[0], _pts[_pts.size() / 2], _pts.back()};
  double angle[3];
  for (int i = 0; i < 3; ++i) {
    c[i] = (c[i] - center).normalized() * radius;
    angle[i] = atan2(c[i][1], c[i][0]);
  }
  if (angle[2] < angle[0])
    angle[2] += PI * 2.;
  if (angle[1] < angle[0])
    angle[1] += PI * 2.;
  if (angle[1] <= angle[2]) { // OK--CCW arc
    double a = angle[2] - angle[0];
    return new Arc(c[0] + center, angle[0] + PI * 0.5, a * radius, 1. / radius);
  } else { // Backwards--CW arc
    for (int i = 0; i < 3; ++i)
      angle[i] = atan2(c[i][1], c[i][0]);
    if (angle[0] <= angle[2])
      angle[0] += PI * 2.;
    if (angle[1] < angle[2])
      angle[1] += PI * 2.;
    assert(angle[1] <= angle[0]);
    double a = angle[0] - angle[2];
    return new Arc(c[0] + center, angle[0] - PI * 0.5, a * radius,
                   -1. / radius);
  }
}

void ClothoidFitter::addPoint(const Vector2d &pt) {
  _pts.push_back(pt);
  if ((int)_pts.size() < 2)
    return;

  const Vector2d &prevPt = _pts[_pts.size() - 2];
  double segmentLength = (pt - prevPt).norm();

  // Update weighted center of mass.
  _centerOfMass = _centerOfMass + (pt + prevPt) * (0.5 * segmentLength);

  // Adjust angle for continuity.
  double angle = atan2(pt[1] - prevPt[1], pt[0] - prevPt[0]);
  if (angle < _prevAngle) // make sure it's not far from the previous angle
    angle += TWOPI * int(0.5 + (_prevAngle - angle) / TWOPI);
  else
    angle -= TWOPI * int(0.5 + (angle - _prevAngle) / TWOPI);
  _prevAngle = angle;

  // Update total curve length.
  double s0 = _totalLength;
  _totalLength += segmentLength;
  double s1 = _totalLength;

  // p(s) = y*s + z
  // y = current angle at s
  // y*s + z = _angIntegral for continuity, or position.
  double z = _angleIntegral - angle * s0;
  _angleIntegral += segmentLength * angle;

  // Accumulate the integral from s0 to s1.
  _rhs += _getRhs(s1, angle, z) - _getRhs(s0, angle, z);
}

ClothoidPtr ClothoidFitter::getCurve() const {
  // The model is a function of the position P(s) = as^3 + bs^2 + cs + d
  // where s is the arc length and a, b, c, d are the coefficients to solve for.
  // We use this function because the derivative of the position is the
  // tangent T(s) = 3as^2 + 2bs + c.  Taking the derivative again gives us the
  // curvature K(s) = 6as + 2b.
  // If we can solve for the coefficients of P(s) to get T(s) then use T'(s)
  // to get K(s)
  // This will approximate the curvature.
  //
  // To solve for P(s), we try to minimize the error:
  //        E(a, b, c, d) = int_0^L (P(s) - P_data(s))^2 ds
  // where P_data(s) is the Position at arc length of the data.
  //
  // We can differentiate the error with respect to a, b, c, d as
  //        dE/da = 2 int_0^L (P(s) - P_data(s)) s^3 ds = 0
  //
  // Rearranging the terms, we get the equation:
  //        2 int_0^L s^3 (as^3 + bs^2 + cs + d) ds = 2 int_0^L s^3 P_data(s) ds
  //        2 int_0^L s^2 (as^3 + bs^2 + cs + d) ds = 2 int_0^L s^2 P_data(s) ds
  //        2 int_0^L s^1 (as^3 + bs^2 + cs + d) ds = 2 int_0^L s^1 P_data(s) ds
  //        2 int_0^L s^0 (as^3 + bs^2 + cs + d) ds = 2 int_0^L s^0 P_data(s) ds
  //
  // These can be set as the matrix equation LHS * [a b c d]^T = RHS

  Matrix4d lhs;

  lhs = _getLhs(_totalLength);
  Vector4d abcd = lhs.inverse() * _rhs;
  return getClothoidWithParams(abcd);
}

ClothoidPtr ClothoidFitter::getCurveWithZeroCurvature(double param) const {
  Matrix<double, 5, 5> lhs;
  Matrix<double, 5, 1> rhs;

  Vector4d constraint;
  constraint << 6 * param, 2, 0, 0;
  // second derivative of ax^3+bx^2+cx+d is 6ax+2b

  // For constrained least squares,
  // lhs is now [A^T A    C]
  //            [ C^T     0]
  lhs << _getLhs(_totalLength), constraint, constraint.transpose(), 0;
  rhs << _rhs, 0;

  Vector4d abcd = (lhs.inverse() * rhs).head<4>();
  return getClothoidWithParams(abcd);
}

ClothoidPtr
ClothoidFitter::getClothoidWithParams(const Eigen::Vector4d &abcd) const {
  // Theta(s) = as^3 + bs^2 + cs + d
  // Theta'(s) = 3as^2 + 2bs + c
  Vector3d abc = Vector3d(abcd[0] * 3, abcd[1] * 2, abcd[2]);

  double startAngle = abc[2];
  double startCurvature = abc[1];
  double endCurvature = 2 * abc[0] * _totalLength + abc[1];

  // now compute the center of mass of the clothoid at 0
  double x, y;

  if (fabs(abc[0]) > 1e-8) // if it's a real clothoid
  {
    // The following comes from the expression that Mathematica generates with:
    // Integrate[
    //  Integrate[ {Cos[a x^2 + b x + c], Sin[a x^2 + b x + c]}, {x, 0, t} ],
    //  {t, 0, s}]/s
    bool negative = false;
    if (abc[0] < 0) {
      negative = true;
      abc = -abc;
    }
    double a = abc[0], b = abc[1], c = abc[2];
    double s = _totalLength;

    double invRtApi2 = 1. / sqrt(0.5 * PI * a);

    double f1s, f2s, f1c, f2c;
    fresnelApprox(b * invRtApi2 * 0.5, &f1s, &f1c);
    fresnelApprox((b + 2 * a * s) * invRtApi2 * 0.5, &f2s, &f2c);

    double disc = b * b / (4 * a) - c;
    double sind = sin(disc), cosd = cos(disc);

    y = (cos(c + s * (b + a * s)) - cos(c)) / (2. * a);
    y -= (b + 2 * a * s) * PI * 0.25 *
         (cosd * (f1s - f2s) + sind * (f2c - f1c)) * invRtApi2 / a;
    if (negative)
      y = -y;

    x = (sin(c) - sin(c + s * (b + a * s))) / (2. * a);
    x += (b + 2 * a * s) * PI * 0.25 *
         (-sind * (f1s - f2s) + cosd * (f2c - f1c)) * invRtApi2 / a;

    x /= s;
    y /= s;
  } else {
    double b = abc[1], c = abc[2];
    double s = _totalLength;

    if (fabs(b) < 1e-8) // line
    {
      x = cos(c) * s * 0.5;
      y = sin(c) * s * 0.5;
    } else // arc
    {
      double c1 = cos(c), s1 = sin(c);
      double c2 = cos(c + b * s), s2 = sin(c + b * s);

      x = (c1 - c2) / (b * b * s) - s1 / b;
      y = (s1 - s2) / (b * b * s) + c1 / b;
    }
  }

  return new Clothoid(Vector2d(_centerOfMass[0] / _totalLength - x,
                               _centerOfMass[1] / _totalLength - y),
                      startAngle, _totalLength, startCurvature, endCurvature);
}

Matrix4d ClothoidFitter::_getLhs(double s) {
  // Precompute powers of s: sp[i] = s^i for i=0..7.
  double sp[8] = {1, s, 0, 0, 0, 0, 0, 0};
  for (int i = 2; i < 8; ++i)
    sp[i] = sp[i - 1] * s;

  Matrix4d out;
  // Build a 4x4 matrix where each element is a scaled power of x.
  // The element at (i,j) uses p = 7 - i - j and represents an integrated
  // moment.
  //
  // 2/7 * x^7  2/6 * x^6  2/5 * x^5  2/4 * x^4
  // 2/6 * x^6  2/5 * x^5  2/4 * x^4  2/3 * x^3
  // 2/5 * x^5  2/4 * x^4  2/3 * x^3  2/2 * x^2
  // 2/4 * x^4  2/3 * x^3  2/2 * x^2  2/1 * x^1
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      int p = 7 - i - j;
      out(i, j) = (2.0 / double(p)) * sp[p];
    }
  return out;
}

Vector4d ClothoidFitter::_getRhs(double s, double y, double z) {
  // Precompute powers of x: xp[i] = x^i for i=0..5.
  double xp[6] = {1, s, 0, 0, 0, 0};
  for (int i = 2; i < 6; ++i)
    xp[i] = xp[i - 1] * s;

  Vector4d out;
  // Build the right-hand side vector from weighted sums of y and z.
  out[0] = (2.0 / 5.0) * y * xp[5] + 0.5 * z * xp[4];
  out[1] = 0.5 * y * xp[4] + (2.0 / 3.0) * z * xp[3];
  out[2] = (2.0 / 3.0) * y * xp[3] + z * xp[2];
  out[3] = y * xp[2] + 2.0 * z * xp[1];
  return out;
}

// ============================================================================
//                            Anticlothoid Fitting
// ============================================================================

AnticlothoidPtr AnticlothoidFitter::getCurve() const { return getCurvev2(); }

AnticlothoidPtr AnticlothoidFitter::getCurvev1() const {
  if (_pts.size() < 3) {
    return new Anticlothoid(Vector2d(0, 0), 0, 1, 1, 1);
  }

  double _error;
  if (!verifyCandidate(_error)) {
    return new Anticlothoid(Vector2d(0, 0), 0, 0, 0, 0);
  }

  Eigen::VectorXd arcLengths = calculateArcLengths();
  Eigen::VectorXd radii = calculateRadiusOfCurvatures();

  // We only use interior points
  assert(arcLengths.size() == _pts.size() - 1);
  assert(arcLengths.size() == radii.size() + 1);

  // get Sn = S(tn) - S(t1) = S(pn) - S(p1)
  // get radius of curvature's rn, r1
  // Say we have [p0, p1, p2, p3, p4], [L0, L1, L2, L3], [r0, r1, r2]
  // p1, L0, r0
  // p3, L2, r2
  double s0 = arcLengths[0], sn_1 = arcLengths[arcLengths.size() - 2];
  double Sn = sn_1 - s0;
  double r1 = radii[0], rn = radii[radii.size() - 1];

  double rnr1 = (rn * rn - r1 * r1);
  double a = rnr1 / (2 * Sn); // a = (rn^2 - r1^2)/(2Sn)
  // double t1 = (2 * Sn * r1) / rnr1; // t1 = 2(Sn r1)/(rn^2 - r1^2)
  // double tn = (2 * Sn * rn) / rnr1; // tn = 2(Sn rn)/(rn^2 - r1^2)

#if 0
  // Debugging::get()->printf("Solve for Params: s0:%f, sn:%f, Sn:%f", s0, sn,
  // Sn); Debugging::get()->printf("                : r1:%f, rn:%f, rnr1:%f",
  // r1, rn,
  //                          rnr1);
  // Debugging::get()->printf("                : a:%f, t1:%f, rn:%f", a, t1,
  // tn);

  // Debugging::get()->printf("Predicted Ts:");
  int count = 20;
  VectorXd ts(count);
  for (double i = 0; i < count; ++i) {
    double t_step = (tn - t1) / (count - 2 - 1);
    ts[i] = t1 + t_step * (i - 1);
    // Debugging::get()->printf("%f", ts[i]);
  }

  double radius = a;
  VectorC<Vector2d> pts(ts.size(), NOT_CIRCULAR);
  for (int i = 0; i < ts.size(); ++i) {
    double t = ts[i];
    double x = a * cos(t) + t * a * sin(t);
    double y = a * sin(t) - t * a * cos(t);
    pts[i] = Vector2d(x, y);
  }

  // Debugging::get()->printf("predicted x(t), y(t)");
  for (auto p : pts) {
    // Debugging::get()->printf("%f %f", p[0], p[1]);
  }

#endif
  Vector2d start = _pts[0];
  // This could be approximated better with t instead of discretized?
  double firstAngle = AngleUtils::toRange(
      atan2(_pts[1][1] - _pts[0][1], _pts[1][0] - _pts[0][0]));
  double nextAngle = AngleUtils::toRange(
      atan2(_pts[2][1] - _pts[1][1], _pts[2][0] - _pts[1][0]));

  // Handle crossing the x axis
  double rawDelta = nextAngle - firstAngle;
  // remainder(x, 2π) returns x wrapped into (−π,π]
  double delta = std::remainder(rawDelta, 2 * M_PI);
  double startAngle = firstAngle - delta * 0.5;

  // double deltaAngle = nextAngle - firstAngle;
  // double startAngle = firstAngle - deltaAngle / 2;

  // Construct the anticlothoid using the parameters:
  // start, startAngle, length (Sn), starting radius (r1),
  // and DRADIUS.
  // Debugging::get()->printf("Anticlothoid: P : (% 4.2f, % 4.2f), startAngle :
  // % "
  //                          "4.2f, length:%4.2f, RoC:%4.2f, a : % 4.2f ",
  //                          start[0], start[1], startAngle, Sn, r1, a);
  double length = arcLengths[arcLengths.size() - 1];

  // Find (s, r^2) for s = 0
  double rnrn = rn * rn;
  double r1r1 = r1 * r1;
  double slope = (rnrn - r1r1) / (sn_1 - s0);
  double dy = slope * s0;
  double r0r0 = r1r1 - dy;
  double r0 = sqrt(abs(r0r0));
  // Debugging::get()->printf("r1: %4.2f, rn: %4.2f, r1r1: %4.2f, rnrn: %4.2f, "
  //                          "slope: %4.2f, dy: %4.2f, r0r0: %4.2f, r0: %4.2f",
  //                          r1, rn, r1r1, rnrn, slope, dy, r0r0, r0);

  // Debugging::get()->printf("start: %4.2f %4.2f, startAngle: %4.2f, length: "
  //                          "%4.2f, r0: %4.2f, a: %4.2f",
  //                          start[0], start[1], startAngle, length, r0, a);
  return new Anticlothoid(start, startAngle, length, r0, a);
}

AnticlothoidPtr AnticlothoidFitter::getCurvev2() const {
  if (_pts.size() < 3) {
    return new Anticlothoid(Vector2d(0, 0), 0, 1, 1, 1);
  }

  double _error;
  if (!verifyCandidate(_error)) {
    return new Anticlothoid(Vector2d(0, 0), 0, 0, 0, 0);
  }

  Eigen::VectorXd arcLengths = calculateArcLengths();
  assert(arcLengths.size() == _pts.size() - 1);

  Eigen::VectorXd radii = calculateRadiusOfCurvatures();
  assert(radii.size() == arcLengths.size() - 1);

  // We only use interior points
  double slope, intercept;
  VectorXd arcLengthshead = arcLengths.head(arcLengths.size() - 1);
  VectorXd radiiSq = radii.array().square();
  _error = fitLineLeastSquares(arcLengthshead, radiiSq, slope, intercept);

  // Sn = S(tn) - S(t1) = S(pn) - S(p1)
  // get radius of curvature's rn, r1
  // Say we have [p1, p2, p3, p4, p5], [s1, s2, s3, s4], [r1, r2, r2]
  // p1 ---- p2 ---- p3 ---- p4 ---- p5
  //     s1      s2      s3      s4
  //         r1      r2      r3
  //
  // p2, L1, r1
  // p4, L3, r3
  double s1 = arcLengths[0], sn_1 = arcLengths[arcLengths.size() - 2];
  double Sn = sn_1 - s1;
  double r1 = sqrt(slope * s1 + intercept);
  double rn = sqrt(slope * sn_1 + intercept);

  double rnr1 = (rn * rn - r1 * r1);
  double a = rnr1 / (2 * Sn); // a = (rn^2 - r1^2)/(2Sn)
  double r0 = sqrt(intercept);
  // double t1 = (2 * Sn * r1) / rnr1; // t1 = 2(Sn r1)/(rn^2 - r1^2)
  // double tn = (2 * Sn * rn) / rnr1; // tn = 2(Sn rn)/(rn^2 - r1^2)
  Debugging::get()->printf("s1: %4.2f, sn_1: %4.2f, Sn: %4.2f, "
                           "r1: %4.2f, rn: %4.2f, rnr1: %4.2f, "
                           "a: %4.2f, r0: %4.2f",
                           s1, sn_1, Sn, r1, rn, rnr1, a, r0);

  Vector2d start = _pts[0];

  double firstAngle = atan2(_pts[1][1] - _pts[0][1], _pts[1][0] - _pts[0][0]);
  double nextAngle = atan2(_pts[2][1] - _pts[1][1], _pts[2][0] - _pts[1][0]);
  firstAngle = AngleUtils::toRange(firstAngle);
  nextAngle = AngleUtils::toRange(nextAngle);
  double rawDelta = nextAngle - firstAngle; // Handle crossing the x axis
  double delta = std::remainder(rawDelta, 2 * M_PI); // wraps into (-π, π]
  double startAngle = firstAngle - delta * 0.5;

  double length = arcLengths[arcLengths.size() - 1];

  // Find (s, r^2) for s = 0
  // double rnrn = rn * rn;
  // double r1r1 = r1 * r1;
  // double slope = (rnrn - r1r1) / (sn_1 - s1);
  // double dy = a * s1;
  // double r0r0 = r1r1 - dy;
  // double r0 = sqrt(abs(r0r0));

  // Debugging::get()->printf("r1: %4.2f, rn: %4.2f, "
  //                          "a: %4.2f, dy: %4.2f, r0r0: %4.2f, r0: %4.2f",
  //                          r1, rn, a, dy, r0r0, r0);

  // Debugging::get()->printf("start: %4.2f %4.2f, startAngle: %4.2f, length: "
  //                          "%4.2f, r0: %4.2f, a: %4.2f",
  //                          start[0], start[1], startAngle, length, r0, a);
  return new Anticlothoid(start, startAngle, length, r0, a);
}

void AnticlothoidFitter::addPoint(const Vector2d &pt) {
  // Similar to ClothoidFitter.
  _pts.push_back(pt);
  if ((int)_pts.size() < 2)
    return;

  const Vector2d &prevPt = _pts[_pts.size() - 2];
  double segmentLength = (pt - prevPt).norm();

  // Update weighted center of mass.
  _centerOfMass = _centerOfMass + (pt + prevPt) * (0.5 * segmentLength);

  // Adjust angle for continuity.
  double angle = atan2(pt[1] - prevPt[1], pt[0] - prevPt[0]);
  if (angle < _prevAngle) // make sure it's not far from the previous angle
    angle += TWOPI * int(0.5 + (_prevAngle - angle) / TWOPI);
  else
    angle -= TWOPI * int(0.5 + (angle - _prevAngle) / TWOPI);
  _prevAngle = angle;

  // Update total curve length.
  double s0 = _totalLength;
  _totalLength += segmentLength;
  double s1 = _totalLength;
}

bool AnticlothoidFitter::verifyCandidate(double &error) const {
  if (_pts.size() < 3) {
    return -1;
  }

  VectorXd arcLengthsRaw = calculateArcLengths();
  VectorXd arcLengths = arcLengthsRaw.head(arcLengthsRaw.size() - 1);
  VectorXd radii = calculateRadiusOfCurvatures();
  assert(arcLengths.size() == radii.size());

  // Perform Least Sqaures fit
  double slope, intercept;
  VectorXd radiiSq = radii.array().square();
  error = fitLineLeastSquares(arcLengths, radiiSq, slope, intercept);
  exportPlotData("./verify.csv", arcLengths, radiiSq, slope, intercept, error);

  // return error < 30; // TODO: Pick error value.
  return true;
}

Eigen::VectorXd AnticlothoidFitter::calculateArcLengths() const {
  // Returns the arc lengths between each point
  if (_pts.size() < 3)
    return Eigen::VectorXd();

  Eigen::VectorXd arcLengths(_pts.size() - 1);
  double accumulated = 0.0;

  for (size_t i = 0; i < _pts.size() - 1; ++i) {
    accumulated += (_pts[i + 1] - _pts[i]).norm();
    arcLengths[i] = accumulated;
  }

  return arcLengths;
}

Eigen::VectorXd AnticlothoidFitter::calculateRadiusOfCurvatures() const {
  Eigen::VectorXd radii;

  if (_pts.size() < 3) {
    return radii; // Not enough points for radius of curvature
  }

  radii = Eigen::VectorXd(_pts.size() - 2);
  for (size_t i = 1; i + 1 < _pts.size(); ++i) {
    Vector2d v1 = _pts[i] - _pts[i - 1];
    Vector2d v2 = _pts[i + 1] - _pts[i];

    double len1 = v1.norm();
    double len2 = v2.norm();
    // Debugging::get()->printf("v1:(%6.2f, %6.2f), v2:(%6.2f, %6.2f)", v1[0],
    //                          v1[1], v2[0], v2[1]);

    if (len1 == 0 || len2 == 0) {
      radii[i - 1] = numeric_limits<double>::infinity(); // Infinite radius =
                                                         // zero curvature
      continue;
    }

    double cos_theta = v1.normalized().dot(v2.normalized());
    cos_theta = clamp(cos_theta, -1.0, 1.0); // Clamp for numerical safety
    double theta = acos(cos_theta);

    double curvature = 2.0 * sin(theta * 0.5) / max(1e-6, sqrt(len1 * len2));

    if (curvature < 1e-10) {
      radii[i - 1] = numeric_limits<double>::infinity(); // Avoid division
                                                         // by near-zero
    } else {
      radii[i - 1] = 1.0 / curvature;
    }
    // Debugging::get()->printf(
    //     "radii[%d]: %6.4f, cos_theta: %6.4f, theta: %6.4f, curvature: %6.4f",
    //     i - 1, radii[i - 1], cos_theta, theta, curvature);
  }

  return radii;
}

double AnticlothoidFitter::fitLineLeastSquares(const VectorXd &x,
                                               const VectorXd &y, double &slope,
                                               double &intercept) const {
  assert(x.size() == y.size());
  const int N = x.size();
  if (N < 2) {
    slope = 0;
    intercept = y.mean(); // degenerate case
    return 0.0;
  }

  // Set up A = [x, 1]
  MatrixXd A(N, 2);
  A.col(0) = x;
  A.col(1) = VectorXd::Ones(N);

  // Solve least squares: A * [slope, intercept] ≈ y
  Vector2d coeffs =
      A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
  slope = coeffs[0];
  intercept = coeffs[1];

  // Compute residuals
  VectorXd prediction = A * coeffs;
  double rmse = std::sqrt((y - prediction).squaredNorm() / N);
  return rmse;
}

void AnticlothoidFitter::exportPlotData(const string &filename,
                                        const Eigen::VectorXd &arcLengths,
                                        const Eigen::VectorXd &values,
                                        const double slope,
                                        const double intercept,
                                        const double error) const {
  ofstream file(filename);
  if (!file.is_open()) {
    cerr << "Could not write to " << filename << endl;
    return;
  }

  file << "s,value,fitted," << slope << "," << intercept << "," << error
       << "\n";
  for (int i = 0; i < arcLengths.size(); ++i) {
    double s = arcLengths[i];
    double v = values[i];
    double fitted = slope * s + intercept;
    file << s << "," << v << "," << fitted << "\n";
  }
  file.close();
}

END_NAMESPACE_Cornu
