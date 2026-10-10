#include "formula_weird/stability/stability.hpp"

#include <algorithm>
#include <cmath>

namespace formula_weird {
namespace {
bool finiteState(const SystemState& state) {
    if (!std::isfinite(state.total_energy)) return false;
    for (const auto& tensor : state.stress_tensors) if (!tensor.allFinite()) return false;
    for (double value : state.plastic_strains) if (!std::isfinite(value) || value < 0.0) return false;
    for (double value : state.temperatures) if (!std::isfinite(value) || value < 0.0) return false;
    for (const auto& value : state.positions) if (!value.allFinite()) return false;
    for (const auto& value : state.velocities) if (!value.allFinite()) return false;
    for (const auto& value : state.accelerations) if (!value.allFinite()) return false;
    for (double mass : state.nodal_masses) if (!std::isfinite(mass) || mass < 0.0) return false;

    // Phase-space arrays must either all be absent (constitutive-only checkpoint)
    // or describe the same number of nodes. This also rejects orphaned velocity,
    // acceleration, or mass arrays when positions is empty.
    const auto node_count = state.positions.size();
    if (state.velocities.size() != node_count ||
        state.accelerations.size() != node_count ||
        state.nodal_masses.size() != node_count) return false;

    for (const auto& shard : state.debris_shards) {
        if (!shard.position.allFinite() || !shard.velocity.allFinite() ||
            !std::isfinite(shard.mass) || shard.mass < 0.0 ||
            !std::isfinite(shard.age) || shard.age < 0.0 ||
            !std::isfinite(shard.lifetime) || shard.lifetime <= 0.0) return false;
    }
    return true;
}
}

StabilityManager::StabilityManager(double energy_growth_tolerance)
    : tolerance_(std::isfinite(energy_growth_tolerance)
                     ? std::max(0.0, energy_growth_tolerance) : 2.0),
      baseline_energy_(0.0), has_snapshot_(false) {}

void StabilityManager::snapshot(const SystemState& state) {
    if (!finiteState(state)) return;
    // Keep the latest accepted checkpoint, not the lowest-energy historical state.
    last_stable_state_ = state;
    baseline_energy_ = state.total_energy;
    has_snapshot_ = true;
}

void StabilityManager::snapshot(const std::vector<Eigen::Matrix3d>& stress,
                                const std::vector<double>& strains,
                                const std::vector<double>& temps,
                                double current_energy) {
    SystemState state;
    state.stress_tensors = stress;
    state.plastic_strains = strains;
    state.temperatures = temps;
    state.total_energy = current_energy;
    snapshot(state);
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
