#ifndef FORMULA_WEIRD_DEBRIS_HPP
#define FORMULA_WEIRD_DEBRIS_HPP

#include <Eigen/Dense>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace formula_weird {

struct DebrisShard {
    Eigen::Vector3d position = Eigen::Vector3d::Zero();
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    double mass = 0.0;
    double age = 0.0;
    double lifetime = 1.0;
};

class DebrisSystem {
public:
    explicit DebrisSystem(double maximum_lifetime = 10.0, std::uint32_t seed = 42);
    bool spawn(const Eigen::Vector3d& origin, const Eigen::Vector3d& initial_velocity,
               std::size_t count, double total_mass);
    void update(double dt, const Eigen::Vector3d& gravity = Eigen::Vector3d(0.0, -9.81, 0.0));
    const std::vector<DebrisShard>& shards() const noexcept;
    void clear() noexcept;
private:
    double maximum_lifetime_;
    std::uint32_t seed_;
    std::vector<DebrisShard> shards_;
};

} // namespace formula_weird
#endif
