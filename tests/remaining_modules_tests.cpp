#include "formula_weird/contact/contact_mapping.hpp"
#include "formula_weird/debris/debris.hpp"
#include "formula_weird/stability/stability.hpp"
#include "formula_weird/time/time_integration.hpp"

#include <cassert>
#include <cmath>
#include <limits>

int main() {
    using namespace formula_weird;

    TimeIntegrator integrator(0.5);
    assert(std::abs(integrator.computeTimeStep(2.0, 4.0) - 0.25) < 1e-12);
    assert(integrator.computeTimeStep(-1.0, 2.0) == 0.0);
    Eigen::Vector3d p(0, 0, 0), v(1, 0, 0), a(0, -1, 0);
    assert(integrator.integrateKinematics(p, v, a, 0.1));
    assert(std::abs(p.y() + 0.01) < 1e-12);
    const auto old_p = p;
    assert(!integrator.integrateKinematics(p, v, a, -1.0));
    assert((p - old_p).norm() == 0.0);

    ContactMapping contact(100.0, 0.5);
    std::vector<Eigen::Vector3d> positions{Eigen::Vector3d(0, -0.1, 0)};
    std::vector<Eigen::Vector3d> velocities{Eigen::Vector3d(1, -1, 0)};
    contact.resolvePlaneContact(positions, velocities, Eigen::Vector3d::Zero(),
                                Eigen::Vector3d::UnitY(), 0.01);
    assert(std::abs(positions[0].y()) < 1e-12);
    assert(positions[0].allFinite() && velocities[0].allFinite());

    StabilityManager stability(2.0);
    SystemState recovered;
    std::vector<Eigen::Matrix3d> stress{Eigen::Matrix3d::Identity()};
    stability.snapshot(stress, {0.1}, {300.0}, 10.0);
    assert(!stability.checkAndRecover(10.0, recovered));
    assert(stability.checkAndRecover(40.0, recovered));
    assert(std::abs(recovered.total_energy - 10.0) < 1e-12);
    assert(stability.checkAndRecover(std::numeric_limits<double>::quiet_NaN(), recovered));

    DebrisSystem debris(4.0, 123);
    assert(!debris.spawn(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 0, 1.0));
    assert(debris.spawn(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 4, 2.0));
    assert(debris.shards().size() == 4);
    debris.update(0.1);
    assert(debris.shards().size() == 4);
    debris.update(-1.0);
    assert(debris.shards().size() == 4);
    debris.clear();
    assert(debris.shards().empty());
}
