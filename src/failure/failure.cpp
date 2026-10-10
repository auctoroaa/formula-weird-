#include "formula_weird/failure/failure.hpp"

#include <Eigen/Eigenvalues>

#include <cmath>
#include <stdexcept>

namespace formula_weird {

FailureModel::FailureModel(double max_plastic_strain,
                           double critical_principal_stress)
    : eps_p_crit_(max_plastic_strain),
      sigma_crit_(critical_principal_stress) {
    if (!std::isfinite(eps_p_crit_) || eps_p_crit_ < 0.0) {
        throw std::invalid_argument(
            "Maximum plastic strain threshold must be finite and non-negative");
    }
    if (!std::isfinite(sigma_crit_) || sigma_crit_ < 0.0) {
        throw std::invalid_argument(
            "Critical principal stress must be finite and non-negative");
    }
}

bool FailureModel::checkFailure(
    double equivalent_plastic_strain,
    const Eigen::Matrix3d& stress_tensor) const {
    if (!std::isfinite(equivalent_plastic_strain) ||
        equivalent_plastic_strain < 0.0) {
        throw std::invalid_argument(
            "Equivalent plastic strain must be finite and non-negative");
    }

    if (!stress_tensor.allFinite()) {
        throw std::invalid_argument("Stress tensor must contain only finite values");
    }

    // A self-adjoint eigensolver assumes symmetry; reject malformed tensors
    // rather than silently using only one triangle of the input.
    constexpr double symmetry_tolerance = 1e-10;
    const double stress_scale = stress_tensor.cwiseAbs().maxCoeff();
    if (!(stress_tensor - stress_tensor.transpose())
             .cwiseAbs()
             .maxCoeff() <= symmetry_tolerance * (1.0 + stress_scale)) {
        throw std::invalid_argument("Stress tensor must be symmetric");
    }

    if (equivalent_plastic_strain >= eps_p_crit_) {
        return true;
    }

    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(stress_tensor);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Failed to compute stress tensor eigenvalues");
    }

    // Eigen returns self-adjoint eigenvalues in ascending order.
    const double max_principal_stress = solver.eigenvalues()(2);
    return max_principal_stress >= sigma_crit_;
}

} // namespace formula_weird
