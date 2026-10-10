#ifndef FORMULA_WEIRD_STABILITY_HPP
#define FORMULA_WEIRD_STABILITY_HPP

#include "formula_weird/debris/debris.hpp"

#include <Eigen/Dense>
#include <vector>

namespace formula_weird {

struct SystemState {
    std::vector<Eigen::Matrix3d> stress_tensors;
    std::vector<double> plastic_strains;
    std::vector<double> temperatures;
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Vector3d> velocities;
    std::vector<Eigen::Vector3d> accelerations;
    std::vector<double> nodal_masses;
    bool element_eroded = false;
    std::vector<DebrisShard> debris_shards;
    double total_energy = 0.0;
};

class StabilityManager {
public:
    // Absolute energy-growth allowance scaled by max(1, |checkpoint energy|).
    explicit StabilityManager(double energy_growth_tolerance = 2.0);

    // Full-state checkpoint. Rejects non-finite values and inconsistent phase-space sizes.
    void snapshot(const SystemState& state);

    // Compatibility overload for constitutive-only callers; kinematic fields remain empty.
    void snapshot(const std::vector<Eigen::Matrix3d>& stress,
                  const std::vector<double>& strains,
                  const std::vector<double>& temps,
                  double current_energy);

    // Returns true and copies the last checkpoint when energy is non-finite or
    // exceeds its checkpoint by the configured allowance.
    bool checkAndRecover(double current_energy, SystemState& state_out);

private:
    double tolerance_;
    double baseline_energy_;
    SystemState last_stable_state_;
    bool has_snapshot_;
};

} // namespace formula_weird
#endif
