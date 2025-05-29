/*--
    Anticlothoid.h

    This file is part of the Cornucopia curve sketching library Fork.
    Copyright (C) 2025 Joonho Kim (joonho@dgp.toronto.edu)

    ==From the previous Copyright==
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

#ifndef CORNUCOPIA_ANTICLOTHOID_H_INCLUDED
#define CORNUCOPIA_ANTICLOTHOID_H_INCLUDED

#include "CurvePrimitive.h"
#include "defs.h"

NAMESPACE_Cornu

    CORNU_SMART_FORW_DECL(Anticlothoid);

class Anticlothoid : public CurvePrimitive {
public:
  enum Param { X = 0, Y, ANGLE, LENGTH, RADIUS, ARADIUS };

  Anticlothoid() {} // uninitialized
  // Note: For a canonical anticlothoid, the endRadius is expected to equal
  // startRadius + LENGTH.
  Anticlothoid(const Vec &start, double startAngle, double length,
               double startRadius, double a);

  // overrides
  void eval(double s, Vec *pos, Vec *der = NULL, Vec *der2 = NULL) const;

  double project(const Vec &point) const;

  double angle(double s) const;
  double curvature(double s) const;
  double radiusOfCurvature(double s) const;

  double startRadius() const { return _params[RADIUS]; }
  double endRadius() const {
    double R0 = _params[RADIUS];
    double a = _params[ARADIUS];
    double t0 = R0 / a;
    double s0 = 0.5 * a * (t0 * t0); // arc-length at t0
    double sn = s0 + _params[LENGTH];
    return sqrt(2 * a * sn);
  }

  PrimitiveType getType() const { return ANTICLOTHOID; }

  void trim(double sFrom, double sTo);
  void flip();
  CurvePrimitivePtr clone() const {
    AnticlothoidPtr out = new Anticlothoid();
    out->setParams(_params);
    return out;
  }
  void derivativeAt(double s, ParamDer &out, ParamDer &outTan) const;
  void derivativeAtEnd(int continuity, EndDer &out) const;

  class _AnticlothoidProjector // internal singleton class
  {
  public:
    virtual double project(const Vec &pt, double from, double to) const = 0;
  };

  double getA() const { return _params[ARADIUS]; }
  bool getSandT(const double s, double &retS, double &retT) const;

  Vec getTransform(Eigen::Matrix2d &mat, Vec &shift) const {
    mat = _mat;
    shift = _startShift;
    return _params.head<2>();
  }

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
protected:
  // override
  void _paramsChanged();
  bool isValidImpl() const;

private:
  bool isFlipped() const { return _params[RADIUS] < 0; }
  double flipCoef() const { return isFlipped() ? -1.0 : 1.0; }

  Vec _startShift; // translation component of transformation from canonical
                   // anticlothoid
  Eigen::Matrix2d
      _mat; // rotation (and possibly reflection) component of transformation

  class _AnticlothoidProjectorImpl;
  static _AnticlothoidProjector *
  _anticlothoidProjector(); // projects onto a generic clothoid
};

END_NAMESPACE_Cornu

#endif // CORNUCOPIA_ANTICLOTHOID_H_INCLUDED
