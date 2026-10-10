#ifndef FORMULA_WEIRD_STABILITY_HPP
#define FORMULA_WEIRD_STABILITY_HPP

#include <Eigen/Dense>
#include <vector>

namespace formula_weird {

struct SystemState {
    std::vector<Eigen::Matrix3d> stress_tensors;
    std::vector<double> plastic_strains;
    std::vector<double> temperatures;
    double total_energy = 0.0;
};

class StabilityManager {
public:
    // tolerance is a non-negative relative growth allowance; a small absolute
    // floor prevents near-zero baselines from causing numerical false alarms.
    explicit StabilityManager(double energy_growth_tolerance = 2.0);

    // Stores the lowest-energy finite state observed (the recovery checkpoint).
    void snapshot(const std::vector<Eigen::Matrix3d>& stress,
                  const std::vector<double>& strains,
                  const std::vector<double>& temps,
                  double current_energy);

    // Returns true and restores state_out when energy is non-finite or exceeds
    // the checkpoint by the configured relative allowance plus absolute floor.
    bool checkAndRecover(double current_energy, SystemState& state_out);

private:
    double tolerance_;
    double baseline_energy_;
    SystemState last_stable_state_;
    bool has_snapshot_;
};

} // namespace formula_weird
#endif
