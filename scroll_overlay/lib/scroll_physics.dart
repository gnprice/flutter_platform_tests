import 'dart:math' as math;

import 'package:flutter/widgets.dart';

class DecicScrollPhysics extends ClampingScrollPhysics {
  const DecicScrollPhysics({super.parent});

  @override
  ClampingScrollPhysics applyTo(ScrollPhysics? ancestor) {
    return DecicScrollPhysics(parent: buildParent(ancestor));
  }

  @override
  Simulation? createBallisticSimulation(ScrollMetrics position, double velocity) {
    final base = super.createBallisticSimulation(position, velocity);
    if (base is! ClampingScrollSimulation)
      return base;
    return DecicScrollSimulation(
      // TODO: hack on velocity clamp
        position: position.pixels, velocity: 8000.0 * velocity.sign, tolerance: tolerance);
  }
}

class DecicScrollSimulation extends Simulation {
  DecicScrollSimulation(
      {super.tolerance, required this.position, required this.velocity}) {
    _duration = _flingDuration(velocity);
    _distance = _flingAverageSpeed(velocity) * _duration;
  }

  final double position;
  final double velocity;

  late double _duration;
  late double _distance;

  static const friction = 0.015;

  // See DECELERATION_RATE.
  static final double _kDecelerationRate = math.log(0.78) / math.log(0.9);

  // See computeDeceleration().
  static double _decelerationForFriction(double friction) {
    return friction * 61774.04968;
  }

  // See getSplineFlingDuration(). Returns a value in seconds.
  double _flingDuration(double velocity) {
    // See mPhysicalCoeff
    final double scaledFriction = friction * _decelerationForFriction(0.84);

    // See getSplineDeceleration().
    final double deceleration = math.log(kInflexion * velocity.abs() / scaledFriction);

    return math.exp(deceleration / (_kDecelerationRate - 1.0));
  }

  // The average speed over the whole fling.  This is the ratio of
  // getSplineFlingDistance to getSplineFlingDuration, in per-second units.
  double _flingAverageSpeed(double velocity) {
    return kInflexion * velocity.abs();
  }

  // See INFLEXION.
  static const kInflexion = 0.35;

  static const kBend = 1 - (3 * kInflexion / 2);

  // Compute x given t, using the approximation.
  static double _flingDistancePenetration(double t) {
    return -41.911 * t*t*t*t*t*t*t*t*t*t + 236.515 * t*t*t*t*t*t*t*t*t + -553.947 * t*t*t*t*t*t*t*t + 686.393 * t*t*t*t*t*t*t + -466.693 * t*t*t*t*t*t + 150.619 * t*t*t*t*t + -12.834 * t*t*t + 2.857 * t;
  }

  // Compute dx/dt given t, using the approximation.
  static double _flingVelocityPenetration(double t) {
    return -419.106 * t*t*t*t*t*t*t*t*t + 2128.634 * t*t*t*t*t*t*t*t + -4431.573 * t*t*t*t*t*t*t + 4804.752 * t*t*t*t*t*t + -2800.157 * t*t*t*t*t + 753.095 * t*t*t*t + -38.502 * t*t + 2.857;
  }

  @override
  double x(double time) {
    final double t = clampDouble(time / _duration, 0.0, 1.0);
    return position + _distance * _flingDistancePenetration(t) * velocity.sign;
  }

  @override
  double dx(double time) {
    final double t = clampDouble(time / _duration, 0.0, 1.0);
    return _distance * _flingVelocityPenetration(t) * velocity.sign / _duration;
  }

  @override
  bool isDone(double time) {
    return time >= _duration;
  }
}


class CubicCubicScrollPhysics extends ClampingScrollPhysics {
  const CubicCubicScrollPhysics({super.parent});

  @override
  ClampingScrollPhysics applyTo(ScrollPhysics? ancestor) {
    return CubicCubicScrollPhysics(parent: buildParent(ancestor));
  }

  @override
  Simulation? createBallisticSimulation(ScrollMetrics position, double velocity) {
    final base = super.createBallisticSimulation(position, velocity);
    if (base is! ClampingScrollSimulation)
      return base;
    return CubicCubicScrollSimulation(
      // TODO: We should clamp the velocity differently.  Currently we clamp
      //   the full 2D velocity to 8000px/s, and then take the component in the
      //   scrollable axis.  This has the effect that if the user moves their
      //   finger diagonally at high speed, we'll end up clamping not just down
      //   to maxFlingVelocity but to a much lower speed.  For example a swipe
      //   at (6240, 8320) px/s interpreted as a vertical scroll will produce a
      //   velocity of not 8000 px/s but 6400 px/s.
      //
      // The relevant code is in [DragGestureRecognizer._checkEnd], calling
      // [Velocity.clampMagnitude].
      //
      // For a quick crude demo here, we just assume every fling is a max-speed
      // fling.  The result: any fast fling lands at pixel-perfect agreement
      // between the Flutter and Android sides!
        position: position.pixels, velocity: 8000.0 * velocity.sign, tolerance: tolerance);
  }
}

class CubicCubicScrollSimulation extends Simulation {
  CubicCubicScrollSimulation(
      {super.tolerance, required this.position, required this.velocity}) {
    _duration = _flingDuration(velocity);
    _distance = _flingAverageSpeed(velocity) * _duration;
  }

  final double position;
  final double velocity;

  late double _duration;
  late double _distance;

  static const friction = 0.015;

  // See DECELERATION_RATE.
  static final double _kDecelerationRate = math.log(0.78) / math.log(0.9);

  // See computeDeceleration().
  static double _decelerationForFriction(double friction) {
    return friction * 61774.04968;
  }

  // See getSplineFlingDuration(). Returns a value in seconds.
  double _flingDuration(double velocity) {
    // See mPhysicalCoeff
    final double scaledFriction = friction * _decelerationForFriction(0.84);

    // See getSplineDeceleration().
    final double deceleration = math.log(kInflexion * velocity.abs() / scaledFriction);

    return math.exp(deceleration / (_kDecelerationRate - 1.0));
  }

  // The average speed over the whole fling.  This is the ratio of
  // getSplineFlingDistance to getSplineFlingDuration, in per-second units.
  double _flingAverageSpeed(double velocity) {
    return kInflexion * velocity.abs();
  }

  // See INFLEXION.
  static const kInflexion = 0.35;

  static const kBend = 1 - (3 * kInflexion / 2);

  // The Android scroll physics has the equation:
  //   t = kBend * z**3 + (1-kBend) * z    (equation 1a)
  //   x =  -1/2 * z**3 +    3/2    * z    (equation 1b)
  // where t is time as a fraction of [_duration] (so 0 <= t <= 1);
  // z is an abstract intermediate variable, with 0 <= z <= 1;
  // and x is the position of the scroll fling, scaled to run from 0 to 1.
  // See the initialization of SPLINE_POSITION (where t, z, x are called
  // alpha, x, and SPLINE_POSITION[i] respectively.)
  //
  // In order to solve for x given t, the Android framework precomputes a table
  // of x at values of t spaced by (1 / NB_SAMPLES) == 1/100 (that table is
  // SPLINE_POSITION), and then does linear interpolation on the values in the
  // table.
  //
  // We instead use an approximation of equation 1a, of the form:
  //   z = a3 * t**3 + a2 * t**2 + a1 * t + a0    (equation 2a)
  // where the coefficients a0, a1, a2, a3 are chosen so that the curve
  // both passes through the desired start and end points (0,0) and (1,1),
  // and has the same slope at those points as equation 1a does.
  // Then we can use the value of z to apply equation 1b directly.
  static const kCoeff3 = (4*kBend*kBend - kBend) / ((1-kBend) * (1 + 2*kBend));
  static const kCoeff2 = -6*kBend*kBend / ((1-kBend) * (1 + 2*kBend));
  static const kCoeff1 = 1 / (1-kBend);

  // Compute z given t, using equation 2a to approximate equation 1a.
  static double _flingProgress(double t) {
    return kCoeff3 * t*t*t + kCoeff2 * t*t + kCoeff1 * t;
  }

  // Compute x given t, using the approximation.
  static double _flingDistancePenetration(double t) {
    final z = _flingProgress(t);
    assert(0 <= z && z <= 1);
    return (- z*z*z + 3*z) / 2;
  }

  // Compute dx/dt given t, using the approximation.
  static double _flingVelocityPenetration(double t) {
    final z = _flingProgress(t);
    assert(0 <= z && z <= 1);
    final dxByDz = (- z*z + 1) * 3/2;
    final dzByDt = 3*kCoeff3 * t*t + 2*kCoeff2 * t + kCoeff1;
    return dxByDz * dzByDt;
  }

  @override
  double x(double time) {
    final double t = clampDouble(time / _duration, 0.0, 1.0);
    return position + _distance * _flingDistancePenetration(t) * velocity.sign;
  }

  @override
  double dx(double time) {
    final double t = clampDouble(time / _duration, 0.0, 1.0);
    return _distance * _flingVelocityPenetration(t) * velocity.sign / _duration;
  }

  @override
  bool isDone(double time) {
    return time >= _duration;
  }
}

// borrowed from framework's foundation/math.dart
double clampDouble(double x, double min, double max) {
  assert(min <= max && !max.isNaN && !min.isNaN);
  if (x < min) {
    return min;
  }
  if (x > max) {
    return max;
  }
  if (x.isNaN) {
    return max;
  }
  return x;
}
