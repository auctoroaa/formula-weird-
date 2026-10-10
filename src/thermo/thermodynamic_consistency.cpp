#include "formula_weird/thermo/thermodynamic_consistency.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace formula_weird {

ThermodynamicConsistency::ThermodynamicConsistency(double density,
                                                   double specific_heat)
    : rho_(density), Cp_(specific_heat) {
    if (!std::isfinite(rho_) || rho_ <= 0.0) {
        throw std::invalid_argument("density must be finite and positive");
    }
    if (!std::isfinite(Cp_) || Cp_ <= 0.0) {
        throw std::invalid_argument("specific heat must be finite and positive");
    }
}

bool ThermodynamicConsistency::verifyDissipation(
    const Eigen::Matrix3d& stress_tensor,
    const Eigen::Matrix3d& plastic_strain_increment,
    double temperature,
    double temperature_increment,
    double dt,
    double& dissipation_out) const {
    dissipation_out = 0.0;

    if (!stress_tensor.allFinite() ||
        !plastic_strain_increment.allFinite() ||
        !std::isfinite(temperature) ||
        !std::isfinite(temperature_increment) ||
        !std::isfinite(dt) || dt <= 0.0 || temperature <= 0.0) {
        return false;
    }

    // This module uses the supplied engineering approximation:
    // D = sigma : (Delta-epsilon_p / dt)
    //     - max(0, rho * Cp * Delta-T / dt).
    // The thermal term is a modeling assumption, not a general identity for
    // Clausius-Duhem dissipation; see the class-level API comment.
    long double plastic_work_increment = 0.0L;
    for (Eigen::Index i = 0; i < 3; ++i) {
        for (Eigen::Index j = 0; j < 3; ++j) {
            plastic_work_increment +=
                static_cast<long double>(stress_tensor(i, j)) *
                static_cast<long double>(plastic_strain_increment(i, j));
        }
    }

    const long double mechanical_rate =
        plastic_work_increment / static_cast<long double>(dt);
    const long double thermal_rate =
        static_cast<long double>(rho_) *
        static_cast<long double>(Cp_) *
        static_cast<long double>(temperature_increment) /
        static_cast<long double>(dt);
    const long double dissipation =
        mechanical_rate - std::max(0.0L, thermal_rate);

    if (!std::isfinite(mechanical_rate) || !std::isfinite(thermal_rate) ||
        !std::isfinite(dissipation) ||
        std::abs(dissipation) >
            static_cast<long double>(std::numeric_limits<double>::max())) {
        dissipation_out = 0.0;
        return false;
    }

    dissipation_out = static_cast<double>(dissipation);
    constexpr double tolerance = 1e-6;
    return dissipation_out >= -tolerance;
}

} // namespace formula_weird
