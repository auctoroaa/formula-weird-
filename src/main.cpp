#include "formula_weird/constitutive/constitutive_model.hpp"
#include "formula_weird/failure/failure.hpp"
#include "formula_weird/rigid/rigid_body_sync.hpp"
#include "formula_weird/thermo/thermodynamic_consistency.hpp"
#include "formula_weird/thermo/thermodynamics.hpp"

#include <Eigen/Dense>

#include <exception>
#include <iostream>
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
        const bool dissipative = consistency.verifyDissipation(
            Eigen::Matrix3d::Zero(), Eigen::Matrix3d::Zero(),
            temperature, 0.0, 1.0e-6, dissipation);
        if (!dissipative) {
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

        std::cout << "formula_weird_engine: available-module smoke tests passed\n";
        std::cout << "NOTE: this repository revision does not yet contain Debris, Stability, "
                     "Time, Contact, umbrella API, or the full explicit-dynamics driver.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "formula_weird_engine failed: " << error.what() << '\n';
        return 1;
    }
}
