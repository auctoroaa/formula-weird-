#include "formula_weird/debris/debris.hpp"

#include <algorithm>
#include <cmath>
#include <random>

namespace formula_weird {
DebrisSystem::DebrisSystem(double maximum_lifetime, std::uint32_t seed)
    : maximum_lifetime_(std::isfinite(maximum_lifetime)
                            ? std::max(0.1, maximum_lifetime) : 10.0),
      seed_(seed) {}

bool DebrisSystem::spawn(const Eigen::Vector3d& origin,
                         const Eigen::Vector3d& initial_velocity,
                         std::size_t count, double total_mass) {
    if (!origin.allFinite() || !initial_velocity.allFinite() || count == 0 ||
        !std::isfinite(total_mass) || total_mass <= 0.0) return false;
    const double shard_mass = total_mass / static_cast<double>(count);
    if (!std::isfinite(shard_mass) || shard_mass <= 0.0) return false;

    std::mt19937 rng(seed_++);
    std::uniform_real_distribution<double> offset(-1.0, 1.0);
    std::uniform_real_distribution<double> life(0.5 * maximum_lifetime_, maximum_lifetime_);
    std::vector<DebrisShard> pending;
    pending.reserve(count);
    for (std::size_t i = 0; i < count; ++i) {
        Eigen::Vector3d direction(offset(rng), offset(rng), offset(rng));
        if (direction.squaredNorm() < 1e-12) direction = Eigen::Vector3d::UnitY();
        direction.normalize();
        DebrisShard shard;
        shard.position = origin;
        shard.velocity = initial_velocity + direction * (0.5 + 2.0 * std::abs(offset(rng)));
        shard.mass = shard_mass;
        shard.lifetime = life(rng);
        pending.push_back(shard);
    }
    shards_.insert(shards_.end(), pending.begin(), pending.end());
    return true;
}

void DebrisSystem::update(double dt, const Eigen::Vector3d& gravity) {
    if (!std::isfinite(dt) || dt <= 0.0 || !gravity.allFinite()) return;
    for (auto& shard : shards_) {
        if (!shard.position.allFinite() || !shard.velocity.allFinite() ||
            !std::isfinite(shard.age) || !std::isfinite(shard.lifetime)) {
            shard.age = maximum_lifetime_;
            continue;
        }
        const Eigen::Vector3d next_velocity = shard.velocity + gravity * dt;
        const Eigen::Vector3d next_position = shard.position + next_velocity * dt;
        const double next_age = shard.age + dt;
        if (!next_velocity.allFinite() || !next_position.allFinite() || !std::isfinite(next_age)) {
            shard.age = maximum_lifetime_;
            continue;
        }
        shard.velocity = next_velocity;
        shard.position = next_position;
        shard.age = next_age;
    }
    shards_.erase(std::remove_if(shards_.begin(), shards_.end(),
        [this](const DebrisShard& shard) { return shard.age >= shard.lifetime; }), shards_.end());
}
const std::vector<DebrisShard>& DebrisSystem::shards() const noexcept { return shards_; }
void DebrisSystem::clear() noexcept { shards_.clear(); }
} // namespace formula_weird
