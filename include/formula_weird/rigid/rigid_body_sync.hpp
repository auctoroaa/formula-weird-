#ifndef FORMULA_WEIRD_RIGID_BODY_SYNC_HPP
#define FORMULA_WEIRD_RIGID_BODY_SYNC_HPP

#include <Eigen/Dense>
#include <vector>

namespace formula_weird {

struct Node {
    Eigen::Vector3d position;
    Eigen::Vector3d velocity;
    double mass;
};

class RigidBodySync {
public:
    RigidBodySync();

    // Recomputes center of mass, linear momentum, and inertia about the center of mass.
    // Nodes with zero mass are ignored; negative or non-finite masses are rejected.
    // If no positive-mass nodes exist, all output values are reset to zero.
    void updateRigidBody(const std::vector<Node>& nodes,
                         Eigen::Vector3d& center_of_mass,
                         Eigen::Matrix3d& inertia_tensor,
                         Eigen::Vector3d& total_momentum);

private:
    Eigen::Vector3d computeCenterOfMass(const std::vector<Node>& nodes,
                                        double total_mass) const;
    Eigen::Matrix3d computeInertiaTensor(const std::vector<Node>& nodes,
                                         const Eigen::Vector3d& com) const;
};

} // namespace formula_weird

#endif // FORMULA_WEIRD_RIGID_BODY_SYNC_HPP
