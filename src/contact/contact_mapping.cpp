#include "formula_weird/contact/contact_mapping.hpp"

#include <algorithm>
#include <cmath>

namespace formula_weird {
namespace {
bool finiteVector(const Eigen::Vector3d& v) { return v.allFinite(); }
}

ContactMapping::ContactMapping(double penalty_stiffness, double friction_coefficient)
    : stiffness_(std::isfinite(penalty_stiffness) ? std::max(0.0, penalty_stiffness) : 0.0),
      friction_(std::isfinite(friction_coefficient)
                    ? std::clamp(friction_coefficient, 0.0, 1.0) : 0.0) {}

void ContactMapping::resolvePlaneContact(std::vector<Eigen::Vector3d>& positions,
                                         std::vector<Eigen::Vector3d>& velocities,
                                         const Eigen::Vector3d& plane_point,
                                         const Eigen::Vector3d& plane_normal,
                                         double dt) {
    if (positions.size() != velocities.size() || !std::isfinite(dt) || dt <= 0.0 ||
        !finiteVector(plane_point) || !finiteVector(plane_normal)) return;

    const double norm = plane_normal.norm();
    if (!std::isfinite(norm) || norm <= 1e-12) return;
    const Eigen::Vector3d n = plane_normal / norm;

    for (std::size_t i = 0; i < positions.size(); ++i) {
        auto& p = positions[i];
        auto& v = velocities[i];
        if (!finiteVector(p) || !finiteVector(v)) continue;

        const double signed_distance = (p - plane_point).dot(n);
        if (!std::isfinite(signed_distance) || signed_distance >= 0.0) continue;

        const double penetration = -signed_distance;
        const Eigen::Vector3d normal_force = stiffness_ * penetration * n;
        const Eigen::Vector3d normal_velocity = v.dot(n) * n;
        const Eigen::Vector3d tangent_velocity = v - normal_velocity;
        const double tangent_speed = tangent_velocity.norm();
        Eigen::Vector3d friction_force = Eigen::Vector3d::Zero();

        // Unit-mass penalty model: Coulomb limit is mu * |F_normal|.
        if (std::isfinite(tangent_speed) && tangent_speed > 1e-12) {
            const double friction_limit = friction_ * normal_force.norm();
            const double braking_force = tangent_speed / dt;
            const double magnitude = std::min(friction_limit, braking_force);
            friction_force = -magnitude * tangent_velocity / tangent_speed;
        }

        const Eigen::Vector3d candidate_velocity = v + (normal_force + friction_force) * dt;
        const Eigen::Vector3d candidate_position = p + penetration * n;
        if (finiteVector(candidate_velocity) && finiteVector(candidate_position)) {
            v = candidate_velocity;
            p = candidate_position;
        }
    }
}
} // namespace formula_weird
