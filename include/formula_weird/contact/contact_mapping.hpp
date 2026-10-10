#ifndef FORMULA_WEIRD_CONTACT_MAPPING_HPP
#define FORMULA_WEIRD_CONTACT_MAPPING_HPP

#include <Eigen/Dense>
#include <vector>

namespace formula_weird {

// Penalty contact for unit-mass particles. Units must be consistent.
// This is a simple plane-contact model, not a general rigid-body solver.
class ContactMapping {
public:
    ContactMapping(double penalty_stiffness, double friction_coefficient);

    void resolvePlaneContact(std::vector<Eigen::Vector3d>& positions,
                             std::vector<Eigen::Vector3d>& velocities,
                             const Eigen::Vector3d& plane_point,
                             const Eigen::Vector3d& plane_normal,
                             double dt);

private:
    double stiffness_;
    double friction_;
};

} // namespace formula_weird
#endif
