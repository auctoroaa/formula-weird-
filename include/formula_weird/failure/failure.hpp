#ifndef FORMULA_WEIRD_FAILURE_FAILURE_HPP
#define FORMULA_WEIRD_FAILURE_FAILURE_HPP

#include <Eigen/Dense>

namespace formula_weird {

class FailureModel {
public:
    // Thresholds must be finite and non-negative. Stress is in the same units
    // as the supplied Cauchy stress tensor; strain is dimensionless.
    FailureModel(double max_plastic_strain, double critical_principal_stress);

    // Returns true if either the accumulated equivalent plastic strain reaches
    // its threshold or the maximum principal stress reaches the Rankine limit.
    // Throws std::invalid_argument for non-finite/negative strain or a
    // non-finite/non-symmetric stress tensor.
    bool checkFailure(double equivalent_plastic_strain,
                      const Eigen::Matrix3d& stress_tensor) const;

private:
    double eps_p_crit_; // Critical accumulated equivalent plastic strain
    double sigma_crit_; // Critical maximum principal tensile stress
};

} // namespace formula_weird

#endif // FORMULA_WEIRD_FAILURE_FAILURE_HPP
