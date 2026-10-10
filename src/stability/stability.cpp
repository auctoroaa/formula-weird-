#include "formula_weird/stability/stability.hpp"

#include <algorithm>
#include <cmath>

namespace formula_weird {
namespace {
bool finiteState(const std::vector<Eigen::Matrix3d>& stress,
                 const std::vector<double>& strains,
                 const std::vector<double>& temps) {
    for (const auto& tensor : stress) if (!tensor.allFinite()) return false;
    for (double value : strains) if (!std::isfinite(value)) return false;
    for (double value : temps) if (!std::isfinite(value)) return false;
    return true;
}
}

StabilityManager::StabilityManager(double energy_growth_tolerance)
    : tolerance_(std::isfinite(energy_growth_tolerance)
                     ? std::max(0.0, energy_growth_tolerance) : 2.0),
      baseline_energy_(0.0), has_snapshot_(false) {}

void StabilityManager::snapshot(const std::vector<Eigen::Matrix3d>& stress,
                                const std::vector<double>& strains,
                                const std::vector<double>& temps,
                                double current_energy) {
    if (!std::isfinite(current_energy) || !finiteState(stress, strains, temps)) return;
    if (has_snapshot_ && current_energy >= baseline_energy_) return;

    last_stable_state_.stress_tensors = stress;
    last_stable_state_.plastic_strains = strains;
    last_stable_state_.temperatures = temps;
    last_stable_state_.total_energy = current_energy;
    baseline_energy_ = current_energy;
    has_snapshot_ = true;
}

bool StabilityManager::checkAndRecover(double current_energy, SystemState& state_out) {
    if (!has_snapshot_) return false;
    const double scale = std::max(1.0, std::abs(baseline_energy_));
    const double allowed_energy = baseline_energy_ + tolerance_ * scale;
    if (!std::isfinite(current_energy) || current_energy > allowed_energy) {
        state_out = last_stable_state_;
        return true;
    }
    return false;
}
} // namespace formula_weird
