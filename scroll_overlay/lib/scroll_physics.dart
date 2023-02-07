import 'dart:math' as math;

import 'package:flutter/widgets.dart';

class PowerLawScrollPhysics extends ClampingScrollPhysics {
  const PowerLawScrollPhysics({super.parent});

  @override
  PowerLawScrollPhysics applyTo(ScrollPhysics? ancestor) {
    return PowerLawScrollPhysics(parent: buildParent(ancestor));
  }

  @override
  Simulation? createBallisticSimulation(ScrollMetrics position, double velocity) {
    final base = super.createBallisticSimulation(position, velocity);
    if (base is! ClampingScrollSimulation)
      return base;
    return PowerLawScrollSimulation(
        position: position.pixels, velocity: velocity, tolerance: tolerance);
  }

  @override
  Simulation? updateBallisticSimulation(Simulation oldSimulation, ScrollMetrics position, double time) {
    return createBallisticSimulation(position, oldSimulation.dx(time));
  }
}

/// An implementation of scroll physics that aligns with Android.
///
/// This travels the same total distance on a given fling as the native Android
/// scroll physics does.  The curve has been adjusted in order to make it
/// ballistic, so that the deceleration at any moment is a function only of the
/// current velocity regardless of how long the simulation has been running.
/// This makes it compatible with Flutter's scroll-physics protocol, where a
/// simulation can be restarted many times using only its current velocity.
class PowerLawScrollSimulation extends Simulation {
  PowerLawScrollSimulation(
      {super.tolerance, required this.position, required this.velocity,
        this.offsetTime = 0.0}) {
    _duration = _flingDuration(velocity);
    // debugPrint('simulation: v0=${velocity.toStringAsFixed(1)}, times ${offsetTime.toStringAsFixed(4)}/${_duration.toStringAsFixed(4)}');
  }

  final double position;
  final double velocity;
  final double offsetTime;

  late double _duration;

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

    // This is scaled so that the total distance traveled is the same as on Android.
    return _kDecelerationRate * kInflexion * math.exp(deceleration / (_kDecelerationRate - 1.0));
  }

  // See INFLEXION.
  static const kInflexion = 0.35;

  double _flingDistancePenetration(double t) {
    return 1.0 - math.pow(1.0 - t, _kDecelerationRate);
  }

  double _flingVelocityPenetration(double t) {
    return math.pow(1.0 - t, _kDecelerationRate - 1.0) as double;
  }

  @override
  double x(double time) {
    final double t = clampDouble((time + offsetTime) / _duration, 0.0, 1.0);
    // debugPrint('x: time ${time.toStringAsFixed(4)}/${_duration.toStringAsFixed(4)} -> ${t.toStringAsFixed(4)}'
    //   + ' -> ${_flingDistancePenetration(t).toStringAsFixed(4)}'
    //   + ' -> ${(position + velocity*_duration/_kDecelerationRate * _flingDistancePenetration(t)).toStringAsFixed(1)}');
    return position + velocity*_duration/_kDecelerationRate * _flingDistancePenetration(t);
  }

  @override
  double dx(double time) {
    final double t = clampDouble((time + offsetTime) / _duration, 0.0, 1.0);
    // debugPrint('dx: time ${time.toStringAsFixed(4)}/${_duration.toStringAsFixed(4)} -> ${t.toStringAsFixed(4)}'
    //     + ' -> ${_flingVelocityPenetration(t).toStringAsFixed(4)}'
    //     + ' -> ${(velocity * _flingVelocityPenetration(t)).toStringAsFixed(1)}');
    return velocity * _flingVelocityPenetration(t);
  }

  @override
  bool isDone(double time) {
    return (time + offsetTime) >= _duration;
  }
}


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
        position: position.pixels, velocity: velocity, tolerance: tolerance);
  }

  @override
  Simulation? updateBallisticSimulation(Simulation oldSimulation, ScrollMetrics position, double time) {
    final base = super.updateBallisticSimulation(oldSimulation, position, time);
    if (base is! ClampingScrollSimulation)
      return base;
    if (oldSimulation is! DecicScrollSimulation) {
      return DecicScrollSimulation(
        position: position.pixels,
        velocity: oldSimulation.dx(time),
        tolerance: tolerance,
      );
    }
    return DecicScrollSimulation(
      position: oldSimulation.position,
      velocity: oldSimulation.velocity,
      tolerance: oldSimulation.tolerance,
      offsetTime: time + oldSimulation.offsetTime,
    );
  }
}

class DecicScrollSimulation extends Simulation {
  DecicScrollSimulation(
      {super.tolerance, required this.position, required this.velocity,
      this.offsetTime = 0.0}) {
    _duration = _flingDuration(velocity);
    _distance = _flingAverageSpeed(velocity) * _duration;
  }

  final double position;
  final double velocity;
  final double offsetTime;

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
    return -41.91062021450815 * t*t*t*t*t*t*t*t*t*t + 236.5149267239282 * t*t*t*t*t*t*t*t*t + -553.9466639943882 * t*t*t*t*t*t*t*t + 686.3932037159542 * t*t*t*t*t*t*t + -466.6927670568513 * t*t*t*t*t*t + 150.6189365754758 * t*t*t*t*t + -12.83415860675335 * t*t*t + 2.857142857142857 * t;
  }

  // Compute dx/dt given t, using the approximation.
  static double _flingVelocityPenetration(double t) {
    return -419.1062021450815 * t*t*t*t*t*t*t*t*t + 2128.634340515354 * t*t*t*t*t*t*t*t + -4431.573311955106 * t*t*t*t*t*t*t + 4804.75242601168 * t*t*t*t*t*t + -2800.156602341108 * t*t*t*t*t + 753.094682877379 * t*t*t*t + -38.50247582026005 * t*t + 2.857142857142857;
  }

  @override
  double x(double time) {
    final double t = clampDouble((time + offsetTime) / _duration, 0.0, 1.0);
    return position + _distance * _flingDistancePenetration(t) * velocity.sign;
  }

  @override
  double dx(double time) {
    final double t = clampDouble((time + offsetTime) / _duration, 0.0, 1.0);
    return _distance * _flingVelocityPenetration(t) * velocity.sign / _duration;
  }

  @override
  bool isDone(double time) {
    return (time + offsetTime) >= _duration;
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
        position: position.pixels, velocity: velocity, tolerance: tolerance);
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
