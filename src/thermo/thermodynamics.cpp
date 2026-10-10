#include "formula_weird/thermo/thermodynamics.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace formula_weird {
namespace {
constexpr double kReferenceTemperature = 293.15;
}

Thermodynamics::Thermodynamics(double taylor_quinney_coeff,
                               double heat_capacity,
                               double thermal_softening_coeff)
    : beta_(taylor_quinney_coeff),
      Cp_(heat_capacity),
      m_(thermal_softening_coeff) {
    if (!std::isfinite(beta_) || beta_ < 0.0 || beta_ > 1.0) {
        throw std::invalid_argument("Taylor-Quinney coefficient must be finite and in [0, 1]");
    }
    if (!std::isfinite(Cp_) || Cp_ <= 0.0) {
        throw std::invalid_argument("Heat capacity must be finite and > 0");
    }
    if (!std::isfinite(m_) || m_ < 0.0) {
        throw std::invalid_argument("Thermal softening coefficient must be finite and >= 0");
    }
}

void Thermodynamics::updateThermalState(double plastic_work_increment,
                                        double density,
                                        double& temperature,
                                        double& current_yield_stress,
                                        double base_yield_stress) const {
    if (!std::isfinite(plastic_work_increment) ||
        !std::isfinite(density) || density <= 0.0 ||
        !std::isfinite(temperature) || temperature < 0.0 ||
        !std::isfinite(base_yield_stress) || base_yield_stress < 0.0) {
        throw std::invalid_argument("Thermal state inputs must be finite and physically valid");
    }

    const double plastic_work = std::max(0.0, plastic_work_increment);
    const double delta_temperature =
        (beta_ * plastic_work) / (density * Cp_);
    const double updated_temperature = temperature + delta_temperature;

    if (!std::isfinite(delta_temperature) || !std::isfinite(updated_temperature)) {
        throw std::overflow_error("Temperature update overflowed");
    }

    const double excess_temperature =
        std::max(0.0, updated_temperature - kReferenceTemperature);
    const double softening_factor = std::max(0.0, 1.0 - m_ * excess_temperature);
    const double updated_yield_stress = base_yield_stress * softening_factor;

    if (!std::isfinite(updated_yield_stress)) {
        throw std::overflow_error("Yield-stress update overflowed");
    }

    temperature = updated_temperature;
    current_yield_stress = updated_yield_stress;
}

} // namespace formula_weird
