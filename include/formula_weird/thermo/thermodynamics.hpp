#pragma once

namespace formula_weird {

/// Converts non-negative plastic work per unit volume (J/m^3) into a
/// temperature increment and updates yield stress using linear thermal
/// softening relative to a reference temperature.
class Thermodynamics {
public:
    /// `taylor_quinney_coeff` is the fraction of plastic work converted to heat
    /// (normally in [0, 1]); `heat_capacity` is J/(kg K); `thermal_softening_coeff`
    /// is 1/K. Throws std::invalid_argument for invalid parameters or state.
    Thermodynamics(double taylor_quinney_coeff,
                   double heat_capacity,
                   double thermal_softening_coeff);

    /// `plastic_work_increment`: J/m^3; `density`: kg/m^3;
    /// `temperature`: K; stresses: Pa. Uses 293.15 K as reference.
    /// Negative plastic work is treated as zero. Throws std::invalid_argument
    /// for non-finite/invalid inputs and std::overflow_error on overflow.
    void updateThermalState(double plastic_work_increment,
                            double density,
                            double& temperature,
                            double& current_yield_stress,
                            double base_yield_stress) const;

private:
    double beta_;
    double Cp_;
    double m_;
};

} // namespace formula_weird
