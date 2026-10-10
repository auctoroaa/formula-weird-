#include "formula_weird/constitutive/constitutive_model.hpp"
#include "formula_weird/contact/contact_mapping.hpp"
#include "formula_weird/debris/debris.hpp"
#include "formula_weird/failure/failure.hpp"
#include "formula_weird/rigid/rigid_body_sync.hpp"
#include "formula_weird/stability/stability.hpp"
#include "formula_weird/thermo/thermodynamic_consistency.hpp"
#include "formula_weird/thermo/thermodynamics.hpp"
#include "formula_weird/time/time_integration.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <vector>

int main() {
    try {
        const auto elastic = formula_weird::constitutive::evaluate(
            {210.0e9, 0.3}, {0.0, 0.0, 0.0, 0.0, 0.0, 0.0});
        if (elastic.strain_energy_density != 0.0) {
            std::cerr << "Constitutive smoke test failed\n";
            return 1;
        }

        formula_weird::Thermodynamics thermo(0.9, 500.0, 1.0e-4);
        double temperature = 293.15;
        double yield_stress = 250.0e6;
        thermo.updateThermalState(0.0, 7800.0, temperature, yield_stress, 250.0e6);

        formula_weird::ThermodynamicConsistency consistency(7800.0, 500.0);
        double dissipation = 0.0;
        if (!consistency.verifyDissipation(Eigen::Matrix3d::Zero(),
                Eigen::Matrix3d::Zero(), temperature, 0.0, 1.0e-6, dissipation)) {
            std::cerr << "Thermodynamic consistency smoke test failed\n";
            return 1;
        }

        formula_weird::FailureModel failure(0.2, 500.0e6);
        if (failure.checkFailure(0.0, Eigen::Matrix3d::Zero())) {
            std::cerr << "Failure-model smoke test failed\n";
            return 1;
        }

        formula_weird::RigidBodySync rigid;
        const std::vector<formula_weird::Node> nodes{
            {Eigen::Vector3d(0.0, 0.0, 0.0), Eigen::Vector3d(1.0, 0.0, 0.0), 1.0},
            {Eigen::Vector3d(2.0, 0.0, 0.0), Eigen::Vector3d(-1.0, 0.0, 0.0), 1.0}
        };
        Eigen::Vector3d center_of_mass;
        Eigen::Matrix3d inertia;
        Eigen::Vector3d momentum;
        rigid.updateRigidBody(nodes, center_of_mass, inertia, momentum);
        if (!center_of_mass.allFinite() || !inertia.allFinite() || !momentum.allFinite()) {
            std::cerr << "Rigid-body smoke test failed\n";
            return 1;
        }

        // Force-driven single-particle kinematics: gravity is force/mass.
        formula_weird::TimeIntegrator integrator(0.5);
        Eigen::Vector3d position(0.0, 1.0, 0.0);
        Eigen::Vector3d velocity(0.0, 0.0, 0.0);
        const double mass = 1.0;
        const Eigen::Vector3d gravity(0.0, -9.81, 0.0);
        const Eigen::Vector3d acceleration = gravity; // F / m, with m = 1 kg.
        const double dt = integrator.computeTimeStep(0.01, 5000.0);
        if (dt <= 0.0 || !integrator.integrateKinematics(position, velocity, acceleration, dt)) {
            std::cerr << "Explicit kinematics smoke test failed\n";
            return 1;
        }

        // For this gravity-only particle, mechanical energy is K + m*g*h.
        const double mechanical_energy =
            0.5 * mass * velocity.squaredNorm() + mass * 9.81 * position.y();
        formula_weird::SystemState checkpoint;
        checkpoint.positions = {position};
        checkpoint.velocities = {velocity};
        checkpoint.accelerations = {acceleration};
        checkpoint.nodal_masses = {mass};
        checkpoint.total_energy = mechanical_energy;
        formula_weird::StabilityManager stability(2.0);
        stability.snapshot(checkpoint);
        formula_weird::SystemState recovered;
        if (!stability.checkAndRecover(std::numeric_limits<double>::quiet_NaN(), recovered) ||
            recovered.positions.size() != 1 || recovered.velocities.size() != 1 ||
            (recovered.positions[0] - position).norm() != 0.0 ||
            (recovered.velocities[0] - velocity).norm() != 0.0) {
            std::cerr << "Full phase-space recovery smoke test failed\n";
            return 1;
        }

        formula_weird::ContactMapping contact(100.0, 0.2);
        std::vector<Eigen::Vector3d> contact_positions{Eigen::Vector3d(0.0, -0.01, 0.0)};
        std::vector<Eigen::Vector3d> contact_velocities{Eigen::Vector3d(0.0, -1.0, 0.0)};
        contact.resolvePlaneContact(contact_positions, contact_velocities,
            Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitY(), dt);
        if (!contact_positions[0].allFinite() || !contact_velocities[0].allFinite()) {
            std::cerr << "Contact smoke test failed\n";
            return 1;
        }

        formula_weird::DebrisSystem debris(4.0, 123);
        if (!debris.spawn(Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero(), 2, 1.0) ||
            debris.shards().size() != 2) {
            std::cerr << "Debris smoke test failed\n";
            return 1;
        }
        debris.update(dt);
        if (debris.shards().size() != 2) {
            std::cerr << "Debris update smoke test failed\n";
            return 1;
        }

        std::cout << "formula_weird_engine: available-module smoke tests passed\n";
        std::cout << "Scope: module integration smoke test, not validation of a fully coupled solver.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "formula_weird_engine failed: " << error.what() << '\n';
        return 1;
    }
}
