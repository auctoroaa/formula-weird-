#include "formula_weird/time/time_integration.hpp"

#include <algorithm>
#include <cmath>

namespace formula_weird {
TimeIntegrator::TimeIntegrator(double safety_factor)
    : safety_factor_(std::isfinite(safety_factor)
                         ? std::clamp(safety_factor, 0.1, 1.0) : 0.5) {}

double TimeIntegrator::computeTimeStep(double characteristic_length, double wave_speed) const {
    if (!std::isfinite(characteristic_length) || !std::isfinite(wave_speed) ||
        characteristic_length <= 0.0 || wave_speed <= 0.0) return 0.0;
    const double dt = safety_factor_ * characteristic_length / wave_speed;
    return std::isfinite(dt) && dt > 0.0 ? dt : 0.0;
}

bool TimeIntegrator::integrateKinematics(Eigen::Vector3d& position,
                                         Eigen::Vector3d& velocity,
                                         const Eigen::Vector3d& acceleration,
                                         double dt) const {
    if (!position.allFinite() || !velocity.allFinite() || !acceleration.allFinite() ||
        !std::isfinite(dt) || dt <= 0.0) return false;
    const Eigen::Vector3d next_velocity = velocity + acceleration * dt;
    const Eigen::Vector3d next_position = position + next_velocity * dt;
    if (!next_velocity.allFinite() || !next_position.allFinite()) return false;
    velocity = next_velocity;
    position = next_position;
    return true;
}
} // namespace formula_weird
