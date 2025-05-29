/*--
    AnticlothoidProjector.cpp

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
#include "Arc.h" // assumed to be a class that represents a circular arc
#include <algorithm>
#include <cmath>
#include <deque>

using namespace std;
using namespace Eigen;
NAMESPACE_Cornu

    // ----------------------------------------------------------------------
    // Helper: _ApproxAntiArc
    //
    // This helper class approximates a short segment of the canonical
    // anti-clothoid by constructing a circular arc through three sample points.
    // The canonical anti-clothoid is defined as:
    //
    //    x(t) = a*(cos t + t*sin t),
    //    y(t) = a*(sin t - t*cos t),
    //
    // so that at t = 0 the point is (a,0) and the curve naturally progresses.
    // ----------------------------------------------------------------------
    class _ApproxAntiArc {
public:
  // 'start' and 'length' are in the canonical parameter t.
  // 'a' is the scale parameter from the anti-clothoid.
  _ApproxAntiArc(double start, double length, double a)
      : _start(start), _length(length), _a(a) {
    // Build three sample points:
    Vector2d p[3];
    for (int i = 0; i < 3; ++i) {
      double t = start + 0.5 * length * i;
      p[i][0] = _a * (cos(t) + t * sin(t));
      p[i][1] = _a * (sin(t) - t * cos(t));
    }
    // Construct a circular arc that interpolates p[0], p[1], and p[2].
    _arc = new Arc(p[0], p[1], p[2]);
  }

  // Test if this approximating arc yields a closer projection to pt.
  // If so, update minDistSq and minT (the best t found so far).
  // 'from' and 'to' delimit the allowable t-range.
  bool test(const Vector2d &pt, double &minDistSq, double &minT, double from,
            double to) const {
    // Quick-reject test: compare the squared difference between the
    // arc's radius and the distance from pt to the arc’s center.
    double dSq = (pt - _arc->center()).squaredNorm();
    double rSq = SQR(_arc->radius());
    if (SQR(dSq - rSq) > 2.0 * minDistSq * (dSq + rSq))
      return false;

    // Compute the projected t on the arc (local to this segment).
    double tLocal = _arc->project(pt);
    // Adjust tLocal into the global canonical parameter:
    tLocal = min(max(tLocal + _start, from), to) - _start;

    double distSq = (pt - _arc->pos(tLocal)).squaredNorm();
    if (distSq >= minDistSq)
      return false;

    minT = _start + tLocal;
    minDistSq = distSq;
    return true;
  }

private:
  ArcPtr _arc;
  double _start;
  double _length;
  double _a;
};

// ----------------------------------------------------------------------
// _AnticlothoidProjectorImpl
//
// This class implements a projector for the canonical anti-clothoid.
// It precomputes a deque of _ApproxAntiArc segments over a range of t
// and then uses them (plus two Newton refinement steps) to find the t
// for which the distance from the canonical anti-clothoid to a given
// point is minimized.
// ----------------------------------------------------------------------
class Anticlothoid::_AnticlothoidProjectorImpl
    : public Anticlothoid::_AnticlothoidProjector {
public:
  typedef Vector2d Vec;

  // Construct the projector with a given scale parameter a.
  _AnticlothoidProjectorImpl(double a) : _arcSpacing(0.1), _a(a) {
    double t;
    // Precompute approximating arcs for t in a symmetric range.
    for (t = 0; t * _arcSpacing < 0.9; t += _arcSpacing) {
      _arcs.push_front(_ApproxAntiArc(-t - _arcSpacing, _arcSpacing, _a));
      _arcs.push_back(_ApproxAntiArc(t, _arcSpacing, _a));
    }
    _maxArcParam = t;
  }

  // Given a query point pt and an allowable range of t [from, to],
  // find the value of t that minimizes the distance from the canonical
  // anti-clothoid to pt.
  double project(const Vec &pt, double from, double to) const {
    double minT, minDistSq;
    Vec startPt, endPt;
    // Evaluate the canonical anti-clothoid at the boundaries.
    eval(from, &startPt);
    eval(to, &endPt);

    minT = from;
    minDistSq = (pt - startPt).squaredNorm();
    double distSq = (pt - endPt).squaredNorm();
    if (distSq < minDistSq) {
      minDistSq = distSq;
      minT = to;
    }

    // Determine indices for the precomputed arcs that overlap [from, to].
    int minArcIdx = (int)floor((_maxArcParam + from) / _arcSpacing);
    if (minArcIdx < 0) {
      minArcIdx = 0;
      // For t values below the precomputed range, approximate on the fly.
      double start = from;
      double stop = min(to, -_maxArcParam);
      int cnt = 0;
      while (start + 1e-8 < stop && ++cnt < 100) {
        double len = min(stop - start, -1.0 / start);
        _ApproxAntiArc(start, len, _a).test(pt, minDistSq, minT, from, to);
        start += len;
      }
    }
    int maxArcIdx = (int)ceil((_maxArcParam + to) / _arcSpacing);
    if (maxArcIdx > (int)_arcs.size()) {
      maxArcIdx = (int)_arcs.size();
      // For t values above the precomputed range, approximate on the fly.
      double start = to;
      double stop = max(_maxArcParam, from);
      int cnt = 0;
      while (start - 1e-8 > stop && ++cnt < 100) {
        double len = min(start - stop, 1.0 / start);
        _ApproxAntiArc(start - len, len, _a)
            .test(pt, minDistSq, minT, from, to);
        start -= len;
      }
    }
    // Test against each precomputed arc in the range.
    for (int i = minArcIdx; i < maxArcIdx; ++i)
      _arcs[i].test(pt, minDistSq, minT, from, to);

    // Refine the result with two Newton iterations.
    minT = projectNewton(minT, pt, from, to);
    minT = projectNewton(minT, pt, from, to);
    return minT;
  }

private:
  // Refine the projection using Newton's method.
  double projectNewton(double guess, const Vec &pt, double from,
                       double to) const {
    Vec p, der, der2;
    eval(guess, &p, &der, &der2);
    double dot = der.dot(pt - p);
    double dotDer = der2.dot(pt - p) - der.squaredNorm();
    if (dotDer >= -1e-30) // if the distance would increase, return guess.
      return guess;
    return max(from, min(to, guess - dot / dotDer));
  }

  // Evaluate the canonical anti-clothoid at parameter t.
  // The canonical evaluation uses:
  //    x(t) = a*(cos t + t*sin t),
  //    y(t) = a*(sin t - t*cos t).
  // Optionally, also compute the first and second derivatives with respect to
  // t.
  void eval(double t, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const {
    if (pos) {
      (*pos)[0] = _a * (cos(t) + t * sin(t));
      (*pos)[1] = _a * (sin(t) - t * cos(t));
    }
    if (der || der2) {
      double ct = cos(t), st = sin(t);
      // First derivative with respect to t:
      //    dx/dt = a*t*cos t,
      //    dy/dt = a*t*sin t.
      if (der)
        *der = Vec(_a * t * cos(t), _a * t * sin(t));
      // Second derivative with respect to t:
      //    d^2x/dt^2 = a*(cos t - t*sin t),
      //    d^2y/dt^2 = a*(sin t + t*cos t).
      if (der2)
        *der2 = _a * Vec(cos(t) - t * sin(t), sin(t) + t * cos(t));
    }
  }

  const double _arcSpacing;    // spacing in the canonical parameter t between
                               // approximations.
  deque<_ApproxAntiArc> _arcs; // precomputed approximating arcs.
  double _maxArcParam;         // maximum |t| used in precomputation.
  double _a; // scale parameter for the canonical anti-clothoid.
};

//
// This function returns a pointer to a (static) anti-clothoid projector
// instance. In a full implementation you might wish to construct a projector
// that depends on the anti-clothoid's parameters; here we assume a fixed scale
// a. (Set 'a' to the desired value; for example, you might retrieve it from an
// Anticlothoid instance.)
//
Anticlothoid::_AnticlothoidProjector *Anticlothoid::_anticlothoidProjector() {
  // For demonstration, we set a fixed value for a.
  static double a = 1.0; // <-- Change as needed!
  static Anticlothoid::_AnticlothoidProjector *projector =
      new _AnticlothoidProjectorImpl(a);
  return projector;
}

END_NAMESPACE_Cornu
