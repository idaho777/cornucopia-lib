/*--
    AnticlothoidVerificationTest.cpp

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
#include "Fitter.h"
#include "Test.h"

#include <Eigen/Dense>
#include <Eigen/Geometry> // for Rotation2D
#include <vector>

#include "Debugging.h"
#include "PrimitiveFitUtils.h"

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace std;
using namespace Eigen;
using namespace Cornu;

class AnticlothoidFitterTest : public TestCase {
  using Vec2 = Vector2d;
  struct FittingTest {
    double a;
    vector<double> ts;
    double noise;
    double validationError;
    double fittingError;
    Vec2 T;
    double theta;

    FittingTest(double a, vector<double> ts, double noise,
                double validationError, double fittingError, const Vector2d &T,
                double theta)
        : a(a), ts(ts), noise(noise), validationError(validationError),
          fittingError(fittingError), T(T), theta(theta) {}

    string getName() {
      std::ostringstream oss;
      oss << std::fixed << std::setprecision(4);
      oss << "FittingTest_";
      oss << "s[" << ts.front() << "_" << ts.back() << "]_";
      oss << "T[" << T[0] << "_" << T[1] << "]_";
      oss << "a[" << a << "]_";
      oss << "noise[" << noise << "]_";
      oss << "Rot[" << theta << "]";

      return oss.str();
    }

    void print() {
      std::ostringstream oss;
      oss << "[";
      for (double t : ts)
        oss << " " << setw(3) << t;
      oss << "]";

      Debugging::get()->printf(
          "FittingTest: a:% 4.2f, noise:% 4.2f, VE:% 4.2f, FE:% 4.2f, T:[% "
          "4.2f,% 4.2f], R:% 4.2f, ts:%s",
          a, noise, validationError, fittingError, T[0], T[1], theta,
          oss.str().c_str());
    }

    AnticlothoidPtr getCurve() {
      double t0 = ts[0];
      double tn = ts[ts.size() - 1];
      double x = a * (cos(t0) + t0 * sin(t0));
      double y = a * (sin(t0) - t0 * cos(t0));

      Vec2 start(x, y);
      double angle = t0;
      double length = 0.5 * a * (tn * tn - t0 * t0);
      double r1 = a * t0;

      start = Rotation2D(theta) * start + T;
      angle += theta;

      return new Anticlothoid(start, angle, length, r1, a);
    }
  };

public:
  std::string name() { return "AnticlothoidFitterTests"; }

  void run() {
    string dirname = "./runs/fits-" + today_YYMMDD() + "/csv/";
    createDirectory(dirname);

    vector<FittingTest> tests = setupPoints();
    for (FittingTest test : tests) {
      // test.print();

      vector<Vec2> pts = sample_ac(test);

      AnticlothoidFitter fitter;
      std::for_each(pts.begin(), pts.end(),
                    [&](const Vec2 p) { fitter.addPoint(p); });
      AnticlothoidPtr ptr = fitter.getCurve();
      // AnticlothoidPtr ptr = test.getCurve();

      exportPlotData(dirname + test.getName() + ".csv", test, pts, ptr);

      // double analytical_length =
      //     0.5 * test.a * (test.ts.back() * test.ts.back());
      // double estimated_length = ptr->length();
      // Debugging::get()->printf("Analytical length: %f, Estimated length: %f",
      //                          analytical_length, estimated_length);

      // Debugging::get()->printf("Exported %s", test.getName().c_str());
      // testValidaion(pts, test);
    }
  }

  void testValidation(vector<Vec2> polyline, FittingTest test) {
    // Run a candidate fitting
    AnticlothoidFitter fitter;
    std::for_each(polyline.begin(), polyline.end(),
                  [&](const Vec2 p) { fitter.addPoint(p); });

    // What is the actual least squares fitting
    double actual_error;
    bool pass = fitter.verifyCandidate(actual_error);
    double expected_error = 40; // TODO: Fix this expected error

    // Check if least squares fitting is less than a value.
    CORNU_ASSERT_LT_MSG(actual_error, expected_error,
                        "Points error is " << actual_error
                                           << " but should be < "
                                           << expected_error);

    // Debugging::get()->printf("Actual error %f", actual_error);
  }

  // ==========================================================================
  //                                Helpers
  // ==========================================================================
  vector<FittingTest> setupPoints() {
    vector<FittingTest> tests;

    vector<vector<double>> ts;
    vector<double> t1;
    for (double s = 0; s < 2; s += 0.2)
      t1.push_back(sqrt(2 * s));
    ts.push_back(t1);

    vector<double> t2;
    for (double s = 2; s < 10; s += 0.2)
      t2.push_back(sqrt(2 * s));
    ts.push_back(t2);

    vector<double> t3;
    for (double s = 5; s < 20; s += 0.2)
      t3.push_back(sqrt(2 * s));
    ts.push_back(t3);

    vector<double> t4;
    for (double s = 9; s < 89; s += 0.2)
      t4.push_back(sqrt(2 * s));
    ts.push_back(t4);

    vector<Vector2d> Ts;
    Ts.push_back(Vector2d(0.0, 0.0));
    Ts.push_back(Vector2d(217.0, 300.0));
    Ts.push_back(Vector2d(400.0, -50.0));
    Ts.push_back(Vector2d(-100.0, -200.0));
    Ts.push_back(Vector2d(-800.0, 300.0));

    vector<double> thetas;
    int thetaN = 23;
    for (int i = 0; i < thetaN; ++i) {
      thetas.push_back(i * 2 * M_PI / thetaN);
    }

    vector<double> as;
    as.push_back(1.0);
    // as.push_back(0.4);
    // as.push_back(76);
    as.push_back(280);
    as.push_back(700);

    vector<double> noise;
    noise.push_back(0.0);
    noise.push_back(0.001);
    noise.push_back(0.005);
    // noise.push_back(1.0);
    // noise.push_back(2.0);

    for (double a : as) {
      for (Vector2d T : Ts) {
        for (double theta : thetas) {
          for (double n : noise) {
            for (vector<double> t : ts) {
              double max_noise = n * n * t.size();
              tests.push_back({a, t, n, max_noise, 0.0, T, theta});
            }
          }
        }
      }
    }

    // // Canonical a=1, 2, 0.5
    // tests.push_back({1, t1, 0.0, 0.0, 0.0, Vec2(0.0, 0.0), 0.0});
    // tests.push_back({2, t1, 0.0, 0.0, 0.0, Vec2(0.0, 0.0), 0.0});
    // tests.push_back({0.5, t1, 0.0, 0.0, 0.0, Vec2(0.0, 0.0), 0.0});
    return tests;
  }

  Vec2 sample_ac(double t, double a) {
    double x = a * cos(t) + t * a * sin(t);
    double y = a * sin(t) - t * a * cos(t);
    return Vec2(x, y);
  }

  vector<Vec2> sample_ac(FittingTest ft) {
    double a = ft.a;
    double noise = ft.noise;
    vector<double> ts = ft.ts;
    Vec2 T = ft.T;
    double theta = ft.theta;

    vector<Vec2> pts;
    pts.reserve(ts.size());
    for (double t : ts) {
      Vec2 pt = sample_ac(t, a);
      pt += a * noise * Vec2::Random(); // Noise
      pts.push_back(Rotation2D(theta) * pt + T);
    }
    return pts;
  }

  void createDirectory(const string dirname) {
    // Delete and create directory
    error_code ec;
    filesystem::remove_all(dirname, ec);
    if (ec) {
      std::cerr << "Warning: failed to remove “" << dirname
                << "”: " << ec.message() << "\n";
    }

    filesystem::create_directories(dirname, ec);
    if (ec) {
      std::cerr << "Could not create directory “" << dirname
                << "”: " << ec.message() << "\n";
      // handle failure…
    }
  }

  void exportPlotData(const string filepath, FittingTest ft, vector<Vec2> pts,
                      AnticlothoidPtr acp) {
    // This outputs a file for plot.py to read.
    ofstream file(filepath);
    if (!file.is_open()) {
      cerr << "Could not write to " << filepath << endl;
      return;
    }

    // == Sampled points ==
    file << "Sampled points" << "\n";
    file << "x,y" << "\n";
    for (Vec2 pt : pts) {
      file << pt[0] << "," << pt[1] << "\n";
    }
    file << "\n";

    // == Anticlothoid fitted ==
    // plot circle
    file << "Fitted Anticlothoid" << "\n";
    Matrix2d mat_f;
    Vec2 T_f;
    acp->getTransform(mat_f, T_f);
    double a_f = acp->getA();
    double cos_theta = mat_f(0, 0);
    double sin_theta = mat_f(1, 0);
    double theta_f = atan2(sin_theta, cos_theta);
    file << "startx,starty,a_radius,theta" << "\n";
    file << T_f[0] << "," << T_f[1] << "," << a_f << "," << theta_f << "\n";
    file << "x,y\n";

    // plot fitted curve segment
    double L_f = acp->length();
    constexpr int N = 20;
    for (int step = 0; step <= N; ++step) {
      double s = L_f * (static_cast<double>(step) / N);
      Vec2 pos;
      acp->eval(s, &pos);
      file << pos[0] << "," << pos[1] << "\n";
    }
    file << "\n";

    // == Anticlothoid Actual==
    file << "Actual Anticlothoid" << "\n";
    // plot circle
    Vec2 T_a = ft.T;
    Matrix2d mat_a;
    double theta_a = ft.theta;
    mat_a << cos(theta_a), -sin(theta_a), sin(theta_a), cos(theta_a);

    double a_a = ft.a;
    file << "startx,starty,a_radius,theta" << "\n";
    file << T_a[0] << "," << T_a[1] << "," << a_a << "," << theta_a << "\n";
    file << "x,y\n";

    // plot fitted curve segment
    double t0 = ft.ts[0];
    double tn = ft.ts[ft.ts.size() - 1];
    for (int step = 0; step <= N; ++step) {
      double t = t0 + (tn - t0) * (static_cast<double>(step) / N);
      double x = a_a * (cos(t) + t * sin(t));
      double y = a_a * (sin(t) - t * cos(t));
      Vec2 pos(x, y);
      pos = mat_a * pos + T_a;
      file << pos[0] << "," << pos[1] << "\n";
    }
    file << "\n";

    // Anticlothoid Actual
    file.close();
  }

  static string getTime() {
    using namespace std::chrono;

    auto now = system_clock::now();
    auto now_time_t = system_clock::to_time_t(now);
    auto now_ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;

    std::tm local_tm = *std::localtime(&now_time_t);

    ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << local_tm.tm_hour << std::setw(2)
        << std::setfill('0') << local_tm.tm_min << std::setw(2)
        << std::setfill('0') << local_tm.tm_sec << std::setw(2)
        << std::setfill('0') << now_ms.count(); // hundredths

    return oss.str(); // e.g., "14230547" → 14:23:05.470
  }

  static string today_YYMMDD() {
    // get current time as time_t
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);

    // convert to local tm structure
    std::tm tm;
#if defined(_MSC_VER)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif

    // format into a string
    std::ostringstream oss;
    oss << std::put_time(&tm, "%y%m%d");
    return oss.str();
  }
};

static AnticlothoidFitterTest test;
