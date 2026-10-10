#include "formula_weird/rigid/rigid_body_sync.hpp"

#include <cmath>
#include <stdexcept>

namespace formula_weird {

RigidBodySync::RigidBodySync() = default;

void RigidBodySync::updateRigidBody(const std::vector<Node>& nodes,
                                    Eigen::Vector3d& center_of_mass,
                                    Eigen::Matrix3d& inertia_tensor,
                                    Eigen::Vector3d& total_momentum) {
    double total_mass = 0.0;
    Eigen::Vector3d momentum_sum = Eigen::Vector3d::Zero();

    for (const auto& node : nodes) {
        if (!std::isfinite(node.mass) || node.mass < 0.0) {
            throw std::invalid_argument("Node mass must be finite and non-negative");
        }
        if (!node.position.allFinite() || !node.velocity.allFinite()) {
            throw std::invalid_argument("Node position and velocity must be finite");
        }
        if (node.mass == 0.0) {
            continue;
        }

        total_mass += node.mass;
        momentum_sum += node.mass * node.velocity;
        if (!std::isfinite(total_mass) || !momentum_sum.allFinite()) {
            throw std::overflow_error("Rigid-body mass or momentum calculation overflowed");
        }
    }

    if (total_mass == 0.0) {
        center_of_mass.setZero();
        inertia_tensor.setZero();
        total_momentum.setZero();
        return;
    }

    center_of_mass = computeCenterOfMass(nodes, total_mass);
    if (!center_of_mass.allFinite()) {
        throw std::overflow_error("Center-of-mass calculation overflowed");
    }

    inertia_tensor = computeInertiaTensor(nodes, center_of_mass);
    if (!inertia_tensor.allFinite()) {
        throw std::overflow_error("Inertia-tensor calculation overflowed");
    }

    total_momentum = momentum_sum;
}

Eigen::Vector3d RigidBodySync::computeCenterOfMass(
    const std::vector<Node>& nodes, double total_mass) const {
    Eigen::Vector3d weighted_pos_sum = Eigen::Vector3d::Zero();

    if (total_mass <= 0.0) {
        return Eigen::Vector3d::Zero();
    }

    for (const auto& node : nodes) {
        if (node.mass == 0.0) {
            continue;
        }
        weighted_pos_sum += node.mass * node.position;
        if (!weighted_pos_sum.allFinite()) {
            throw std::overflow_error("Weighted position sum overflowed");
        }
    }

    return weighted_pos_sum / total_mass;
}

Eigen::Matrix3d RigidBodySync::computeInertiaTensor(
    const std::vector<Node>& nodes, const Eigen::Vector3d& com) const {
    Eigen::Matrix3d inertia = Eigen::Matrix3d::Zero();

    for (const auto& node : nodes) {
        if (node.mass == 0.0) {
            continue;
        }

        const Eigen::Vector3d r = node.position - com;
        const double r_squared = r.squaredNorm();
        if (!std::isfinite(r_squared)) {
            throw std::overflow_error("Squared radius overflowed during inertia calculation");
        }

        inertia.noalias() += node.mass *
            (r_squared * Eigen::Matrix3d::Identity() - r * r.transpose());

        if (!inertia.allFinite()) {
            throw std::overflow_error("Inertia-tensor calculation overflowed");
        }
    }

    return inertia;
}

} // namespace formula_weird
