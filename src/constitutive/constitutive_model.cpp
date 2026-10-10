#include "formula_weird/constitutive/constitutive_model.hpp"

#include <cmath>
#include <string>

namespace formula_weird::constitutive {
namespace {

void validate_material(const MaterialParameters& material) {
    if (!std::isfinite(material.youngs_modulus) ||
        material.youngs_modulus <= 0.0) {
        throw std::invalid_argument("Young's modulus must be finite and > 0 Pa");
    }

    if (!std::isfinite(material.poissons_ratio) ||
        material.poissons_ratio <= -1.0 ||
        material.poissons_ratio >= 0.5) {
        throw std::invalid_argument("Poisson's ratio must be finite and in (-1, 0.5)");
    }
}

} // namespace

Response evaluate(const MaterialParameters& material, const Strain6& strain) {
    validate_material(material);

    for (double component : strain) {
        if (!std::isfinite(component)) {
            throw std::invalid_argument("Strain components must all be finite");
        }
    }

    const double E = material.youngs_modulus;
    const double nu = material.poissons_ratio;
    const double mu = E / (2.0 * (1.0 + nu));
    const double lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));

    const double trace = strain[0] + strain[1] + strain[2];

    Response response{};
    response.stress[0] = lambda * trace + 2.0 * mu * strain[0];
    response.stress[1] = lambda * trace + 2.0 * mu * strain[1];
    response.stress[2] = lambda * trace + 2.0 * mu * strain[2];

    // Engineering shear strain convention: tau = mu * gamma.
    response.stress[3] = mu * strain[3]; // yz
    response.stress[4] = mu * strain[4]; // xz
    response.stress[5] = mu * strain[5]; // xy

    // Energy density = 1/2 sigma:epsilon. In engineering-Voigt notation,
    // the shear contribution is tau * gamma (not 2 * tau * gamma).
    double work = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        work += response.stress[i] * strain[i];
    }
    for (std::size_t i = 3; i < 6; ++i) {
        work += response.stress[i] * strain[i];
    }
    response.strain_energy_density = 0.5 * work;

    if (!std::isfinite(response.strain_energy_density)) {
        throw std::overflow_error("Constitutive response overflowed");
    }
    for (double component : response.stress) {
        if (!std::isfinite(component)) {
            throw std::overflow_error("Constitutive stress overflowed");
        }
    }

    return response;
}

} // namespace formula_weird::constitutive
