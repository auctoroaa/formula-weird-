#ifndef FORMULA_WEIRD_THERMODYNAMIC_CONSISTENCY_HPP
#define FORMULA_WEIRD_THERMODYNAMIC_CONSISTENCY_HPP

#include <Eigen/Dense>

namespace formula_weird {

class ThermodynamicConsistency {
public:
    ThermodynamicConsistency(double density, double specific_heat);

    // Checks the implemented local dissipation-rate approximation.
    // This is a diagnostic, not a complete Clausius-Duhem proof: a full
    // check also requires a constitutive free-energy and entropy model.
    bool verifyDissipation(const Eigen::Matrix3d& stress_tensor,
                           const Eigen::Matrix3d& plastic_strain_increment,
                           double temperature,
                           double temperature_increment,
                           double dt,
                           double& dissipation_out) const;

private:
    double rho_;
    double Cp_;
};

} // namespace formula_weird

#endif // FORMULA_WEIRD_THERMODYNAMIC_CONSISTENCY_HPP
