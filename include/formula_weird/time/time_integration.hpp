#ifndef FORMULA_WEIRD_TIME_INTEGRATION_HPP
#define FORMULA_WEIRD_TIME_INTEGRATION_HPP

#include <Eigen/Dense>

namespace formula_weird {

class TimeIntegrator {
public:
    explicit TimeIntegrator(double safety_factor = 0.5);
    double computeTimeStep(double characteristic_length, double wave_speed) const;
    // Semi-implicit Euler; returns false and leaves state unchanged on invalid input.
    bool integrateKinematics(Eigen::Vector3d& position, Eigen::Vector3d& velocity,
                             const Eigen::Vector3d& acceleration, double dt) const;
private:
    double safety_factor_;
};

} // namespace formula_weird
#endif
