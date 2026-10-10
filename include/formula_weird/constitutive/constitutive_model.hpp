#pragma once

#include <array>
#include <stdexcept>

namespace formula_weird::constitutive {

/// Small-strain, isotropic, linear-elastic material parameters (SI units).
struct MaterialParameters {
    double youngs_modulus;       // Pa
    double poissons_ratio;      // dimensionless
};

/// Voigt order: [xx, yy, zz, yz, xz, xy].
/// Shear strain components use engineering shear strain (gamma = 2 * epsilon).
using Strain6 = std::array<double, 6>;
using Stress6 = std::array<double, 6>;

struct Response {
    Stress6 stress{};            // Pa
    double strain_energy_density{}; // J/m^3
};

/// Computes sigma = C : epsilon for isotropic linear elasticity.
/// Throws std::invalid_argument for non-physical material parameters or
/// non-finite strain values.
Response evaluate(const MaterialParameters& material, const Strain6& strain);

} // namespace formula_weird::constitutive
