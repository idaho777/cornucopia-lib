/*--
    Anticlothoid.cpp

    This file is part of the Cornucopia curve sketching library Fork.
    Copyright (C) 2025 Joonho Kim (joonho@dgp.toronto.edu)

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
#include "AngleUtils.h"
#include "Debugging.h"
#include "Eigen/Dense"
#include <cmath>
#include <limits>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

Anticlothoid::Anticlothoid(const Vec &start, double startAngle, double length,
                           double startRadius, double a) {
  // We now store 6 parameters: X, Y, ANGLE, LENGTH, RADIUS, and A.
  // (Previously, the sixth parameter was DRADIUS, but now we use A.)
  // [X], [Y], [ANGLE] are world-space quantities

  _params.resize(numParams());

  _params.head<2>() = start;
  _params[ANGLE] = AngleUtils::toRange(startAngle);
  _params[LENGTH] = length;
  _params[RADIUS] = startRadius;
  _params[ARADIUS] = a; // index 5 now holds the intrinsic parameter a.

  _paramsChanged();
}

void Anticlothoid::_paramsChanged() {
  double L = _params[LENGTH];  // total arc-length of the segment
  double R0 = _params[RADIUS]; // starting radius (r₁)
  double a = _params[ARADIUS]; // our stored intrinsic parameter
  double t0 = R0 / a;          // starting t-value

#define PARAMSDEBUG 0
#if PARAMSDEBUG
  Debugging::get()->printf("Anticlothoid::paramsChanged: L %f R0 %f a %f t0 %f",
                           L, R0, a, t0);
#endif

  // Evaluate the canonical anticlothoid at t0:
  double cosT = cos(t0);
  double sinT = sin(t0);
  Vector2d Q0;
  Q0[0] = a * (cosT + t0 * sinT);
  Q0[1] = a * flipCoef() * (sinT - t0 * cosT);
#if PARAMSDEBUG
  Debugging::get()->printf("Anticlothoid::paramsChanged: Q0 %f %f", Q0[0],
                           Q0[1]);
#endif

  // Set up rotation: we want the world-space tangent at s=0 to be
  // _params[ANGLE]. In the canonical system, the tangent at t0 is (cos t0, sin
  // t0). So, choose a rotation by ( _params[ANGLE] - t0 ).

  double theta = _params[ANGLE];
  if (isFlipped()) {
    theta = AngleUtils::toRange(theta + M_PI);
  }
#if PARAMSDEBUG
  Debugging::get()->printf("Anticlothoid::paramsChanged: theta %f, t0 %f",
                           theta, t0);
#endif
  double rotAngle = theta - (flipCoef() * t0); // negative t is not the angle.
  double cosA = cos(rotAngle);
  double sinA = sin(rotAngle);
  _mat << cosA, -sinA, sinA, cosA;
  // _mat << 1, 0, 0, 1;

  // _params[X] and _params[Y] are the world coordinates of the starting point.
  Vector2d startPos(_params[X], _params[Y]);
  // We want: _startShift + _mat * Q0 = startPos.
#if PARAMSDEBUG
  Debugging::get()->printf("startPos %f %f", startPos[0], startPos[1]);
  Debugging::get()->printf("Q0 %f %f", Q0[0], Q0[1]);
#endif
  _startShift = startPos - _mat * Q0;
}

void Anticlothoid::eval(double s, Vec *pos, Vec *der, Vec *der2) const {
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];

  double trueS, t;
  bool negS = getSandT(s, trueS, t);

  // Also compute Q(t0):
  double ct = cos(t);
  double st = sin(t);
  double Qx = a * (ct + t * st);
  double Qy = a * flipCoef() * (st - t * ct);
  Vector2d canonicalPt(Qx, Qy);

  // Transform to world space:
  if (pos) {
    Vector2d p = _startShift + _mat * canonicalPt;
    (*pos)[0] = p[0];
    (*pos)[1] = p[1];
  }

  // --- First Derivative ---
  // We have t = sqrt(t0^2 + 2s/_a), so: dt/ds = 1 / (_a * t)
  double dt_ds = 1.0 / (a * t);
  // Derivative of Q(t) with respect to t is:
  //   dQ/dt = (a*t*cos t, a*t*sin t)
  double dx_dt = a * t * ct;
  double dy_dt = a * t * st;
  // Thus, derivative of Q(t) with respect to s:
  double dx_ds = dx_dt * dt_ds;
  double dy_ds = dy_dt * dt_ds;
  // Since Q(t0) is constant, the derivative of (Q(t) - Q(t0)) is the same.
  // Therefore, the canonical derivative of P_c(s) is:
  Vector2d canonicalDer(dx_ds, flipCoef() * dy_ds);
  if (der) {
    Vector2d tangent = _mat * canonicalDer;
    (*der)[0] = tangent[0];
    (*der)[1] = tangent[1];
  }

  // --- Second Derivative ---
  // We have dt/ds = 1/(a*t). Differentiate with respect to s:
  //   d^2t/ds^2 = d/ds [1/(a*t)] = -1/(a*t^2) * (dt/ds) = -1/(a*t^2) *
  //   (1/(a*t)) = -1/(a^2 * t^3).
  double d2t_ds2 = -1.0 / (a * a * t * t * t);
  // Second derivative of Q(t) with respect to t is:
  //   d^2Q/dt^2 = (a*(cos t - t*sin t), a*(sin t + t*cos t))
  double d2Qx_dt2 = a * (ct - t * st);
  double d2Qy_dt2 = a * flipCoef() * (st + t * ct);
  // Chain rule: d^2Q/ds^2 = d^2Q/dt^2 * (dt/ds)^2 + dQ/dt * d^2t/ds^2.
  double d2x_ds2 = d2Qx_dt2 * (dt_ds * dt_ds) + dx_dt * d2t_ds2;
  double d2y_ds2 = d2Qy_dt2 * (dt_ds * dt_ds) + dy_dt * d2t_ds2;
  Vector2d canonicalDer2(d2x_ds2, d2y_ds2);
  if (der2) {
    Vector2d secondDer = _mat * canonicalDer2;
    (*der2)[0] = secondDer[0];
    (*der2)[1] = flipCoef() * secondDer[1];
  }
}

double Anticlothoid::angle(double s) const {
#define RADIUSDEBUG 0
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];
  double t0 = R0 / a;
#if RADIUSDEBUG
  Debugging::get()->printf("Angle()::paramsChanged: s %f R0 %f a %f t0 %f", s,
                           R0, a, t0);
#endif
  double trueS, t;
  bool negS = getSandT(s, trueS, t);
#if RADIUSDEBUG
  Debugging::get()->printf("\ttrueS %f, t %f", trueS, t);
  Debugging::get()->printf("\tangle: %f %f", _params[ANGLE], t - t0);
#endif
  if (isFlipped()) {
    t = -t - M_PI;
  }
  // return AngleUtils::toRange(_params[ANGLE] + (t - t0));
  // Get the matrix rotation angle:
  double matAngle = atan2(_mat(1, 0), _mat(0, 0));
#if RADIUSDEBUG
  Debugging::get()->printf("\tt: %f, matAngle: %f", t, matAngle);
#endif
  return AngleUtils::toRange(t + matAngle);
}

double Anticlothoid::curvature(double s) const {
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];
  double t0 = R0 / a;
  double trueS, t;
  bool negS = getSandT(s, trueS, t);
  if (s < 1e-8)
    return 1.0 / R0; // limiting value at s=0
  return 1.0 / (a * t);
}

double Anticlothoid::radiusOfCurvature(double s) const {
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];
  double trueS, t;
  bool negS = getSandT(s, trueS, t);
  return a * t;
}

double Anticlothoid::project(const Vec &point) const {
  // A simple projection by uniformly sampling s in [0, L]
  double bestS = 0;
  double bestDist = std::numeric_limits<double>::infinity();
  int steps = 100;
  double L = _params[LENGTH]; // arc-length of this anticlothoid segment
  for (int i = 0; i <= steps; i++) {
    double s = L * i / steps;
    Vec pos;
    eval(s, &pos, nullptr, nullptr);
    double dx = pos[0] - point[0];
    double dy = pos[1] - point[1];
    double dist = dx * dx + dy * dy;
    if (dist < bestDist) {
      bestDist = dist;
      bestS = s;
    }
  }
  return bestS;
}

void Anticlothoid::trim(double sFrom, double sTo) {
  // We want to trim the segment so that the new segment corresponds
  // to the arc from s = sFrom to s = sTo (with s measured from the current
  // start).
  Vec newStart;
  eval(sFrom, &newStart, nullptr, nullptr);
  // Update the starting tangent angle to be the world tangent at s = sFrom.
  _params[ANGLE] = angle(sFrom);
  _params[RADIUS] = radiusOfCurvature(sFrom);
  // Update LENGTH to the new segment length.
  _params[LENGTH] = sTo - sFrom;
  // Update the starting point.
  _params[X] = newStart[0];
  _params[Y] = newStart[1];
  // We keep DRADIUS unchanged (the per‐unit change in radius remains the same).
  _paramsChanged();
}

void Anticlothoid::flip() {
  // Flip the curve

  // The new start will be the old end.
  Vec newStart;
  eval(_params[LENGTH], &newStart, nullptr, nullptr);

  // New start angle will be the old end angle + 180 degrees.
  double endAngle = angle(_params[LENGTH]);

  // New radius will be the negative old end radius.
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];
  double t0 = R0 / a;
  double truS, t;
  getSandT(_params[LENGTH], truS, t);

  double endRadius = sqrt(2 * a * truS);
#define FLIPDEBUG 0
#if FLIPDEBUG
  Debugging::get()->printf("Anticlothoid::flip: new_start %f %f", newStart[0],
                           newStart[1]);
  Debugging::get()->printf("Anticlothoid::flip: endAngle %f", endAngle);
  Debugging::get()->printf("Anticlothoid::flip: endRadius %f", endRadius);
#endif

  _params[X] = newStart[0];
  _params[Y] = newStart[1];
  _params[ANGLE] = AngleUtils::toRange(endAngle - M_PI);
  _params[RADIUS] = -flipCoef() * endRadius;

  // Length and A's radius stay the same.
  _paramsChanged();
}

// derivativeAt computes the sensitivity of P(s) with respect to each parameter
// _params indices: 0: X, 1: Y, 2: ANGLE, 3: LENGTH, 4: RADIUS, 5: ARADIUS.
void Anticlothoid::derivativeAt(double s, ParamDer &out,
                                ParamDer &outTan) const {
  const int n = numParams(); // expected to be 6
  // Initialize all sensitivities to zero.
  out = ParamDer::Zero(2, n);
  outTan = ParamDer::Zero(2, n);

  // ---------------------------------
  // (1) Sensitivity with respect to X and Y.
  out(0, X) = 1.0; // dP/dX = (1,0)
  out(1, Y) = 1.0; // dP/dY = (0,1)

  // ---------------------------------
  // (2) Evaluate the curve.
  // Compute: t₀ = R₀ / a,  t = sqrt(t₀² + 2s/a)
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];
  double t0 = R0 / a;
  double t = sqrt(t0 * t0 + 2.0 * s / a);

  // Compute Q(t) and Q(t₀)
  double cos_t = cos(t);
  double sin_t = sin(t);
  double cos_t0 = cos(t0);
  double sin_t0 = sin(t0);
  // Q(t) = [ a*(cos t + t*sin t),  a*(sin t - t*cos t) ]
  Vector2d Q_t;
  Q_t[0] = a * (cos_t + t * sin_t);
  Q_t[1] = a * (sin_t - t * cos_t);
  // Q(t₀) = [ a*(cos t₀ + t₀*sin t₀),  a*(sin t₀ - t₀*cos t₀) ]
  Vector2d Q_t0;
  Q_t0[0] = a * (cos_t0 + t0 * sin_t0);
  Q_t0[1] = a * (sin_t0 - t0 * cos_t0);

  // Let the canonical offset be: ΔQ = Q(t) - Q(t₀)
  Vector2d dQ = Q_t - Q_t0;

  // (3) Sensitivity with respect to ANGLE.
  // The world-space rotation matrix M = R(ANGLE - t₀) and its derivative:
  // dM/dANGLE = M * J, with J = [0, -1; 1, 0].
  Matrix2d J;
  J << 0, -1, 1, 0;
  // Then: ∂P/∂ANGLE = M * J * (Q(t) - Q(t₀))
  Vector2d dP_dANGLE = _mat * (J * dQ);
  out.col(ANGLE) = dP_dANGLE;

  // Also, compute the world-space tangent.
  // Since dP/ds = M*[cos t, sin t]ᵀ, we use that for a LENGTH hint.
  Vector2d tangent = _mat * Vector2d(std::cos(t), std::sin(t));
  outTan.col(LENGTH) =
      tangent; // sensitivity hint: change in LENGTH moves along tangent.

  // Additionally, we can set the ANGLE sensitivity for the tangent.
  // A small change in ANGLE rotates the tangent by 90°:
  outTan.col(ANGLE) = Vector2d(-tangent[1], tangent[0]);

  // (4) For simplicity in this minimal example, we leave LENGTH, RADIUS, and
  // DRADIUS sensitivities as zero. In a full implementation, you would compute:
  //   ∂P/∂p = M * ( (∂Q(t)/∂p - ∂Q(t₀)/∂p) ) - M*J*(∂t₀/∂p)*(Q(t)-Q(t₀))
  // for p ∈ {LENGTH, RADIUS, DRADIUS} by differentiating
  //   a = ((R₀ + DRADIUS*LENGTH)² - R₀²) / (2*LENGTH),
  //   t₀ = R₀ / a, and t = sqrt(t₀² + 2s/a).
  // For now, we set these columns to zero.
  // out.col(LENGTH), out.col(RADIUS), out.col(DRADIUS) remain zero.

  // Similarly, for the tangent sensitivities for these parameters, we leave
  // them zero.
}

// derivativeAtEnd computes the sensitivity at s = LENGTH (the end of the
// segment) and adds extra rows for continuity constraints.
void Anticlothoid::derivativeAtEnd(int continuity, EndDer &out) const {
  ParamDer pDer, dummy;
  // Compute the sensitivity at s = LENGTH.
  derivativeAt(_params[LENGTH], pDer, dummy);
  int n = numParams(); // expected to be 6
  // The output EndDer has (2 + continuity) rows and n columns.
  out = EndDer::Zero(2 + continuity, n);
  // Set the first two rows to be the sensitivity at the end.
  out.block(0, 0, 2, n) = pDer;
  // Also, assign the world–space tangent at s = LENGTH to the LENGTH column.
  Vec endTan = der(_params[LENGTH]);
  out.col(LENGTH).head<2>() = endTan;

  // For G1 continuity, we add an extra row.
  if (continuity >= 1) {
    // Here we enforce that a unit change in ANGLE produces a unit change in the
    // tangent direction.
    out(2, ANGLE) = 1.0;
  }
  // For G2 continuity, additional (dummy) values are provided.
  if (continuity >= 2) {
    out(3, LENGTH) = 1.0;
  }
}

bool Anticlothoid::isValidImpl() const {
  if (_params[LENGTH] < 0.)
    return false;
  if (_params[RADIUS] < 0.)
    return false;
  // DRADIUS can be negative (for a decreasing radius), so no check here.
  return true;
}

bool Anticlothoid::getSandT(const double s, double &retS, double &retT) const {
  double R0 = _params[RADIUS];
  double a = _params[ARADIUS];

  double t0 = R0 / a;
  double s0 = 0.5 * a * (t0 * t0); // arc-length at t0
  double trueS = isFlipped() ? s0 - s : s + s0;

  // Check if we go from large R0 to small Rn
  // This means, we're working with negative t values.
  bool negS = trueS < 0;
  if (negS) {
    trueS *= -1;
  }

  double t = flipCoef() * sqrt(2.0 * trueS / a);
  if (negS) {
    t *= -1;
  }
  retS = trueS;
  retT = t;

  return negS;
}

END_NAMESPACE_Cornu
