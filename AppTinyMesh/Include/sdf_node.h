#pragma once

#include <memory>
#include <stdexcept>
#include <utility>

#include "implicits.h"
#include "box.h"
#include <limits>
#include <cmath>
#include <optional>

namespace sdf
{

// Convention: negative inside the object, zero on its surface, positive outside.
class Node : public AnalyticScalarField
{
public:
  virtual ~Node() = default;

  double Value(const Vector& point) const override = 0;
  // Present only when the node can certify a finite box and a distance
  // lower bound outside that box. An arbitrary 1-Lipschitz field is not enough.
  virtual std::optional<::Box> ConservativeBounds() const { return std::nullopt; }
  virtual double LowerBound(const Vector&) const
  {
    return -std::numeric_limits<double>::infinity();
  }
};

using NodePtr = std::shared_ptr<const Node>;

// Associates a subtree with a certified box-distance lower bound.
class BoundedNode final : public Node
{
public:
  BoundedNode(NodePtr child, const ::Box& bounds)
    : child_(std::move(child)), bounds_(bounds)
  {
    if (!child_) throw std::invalid_argument("A bounded SDF node needs a child");
    for (int axis = 0; axis < 3; ++axis)
      if (!std::isfinite(bounds_[0][axis]) || !std::isfinite(bounds_[1][axis]) ||
          bounds_[0][axis] > bounds_[1][axis])
        throw std::invalid_argument("SDF bounds must be finite and ordered");
    const auto certified = child_->ConservativeBounds();
    if (!certified)
      throw std::invalid_argument("This SDF subtree has no certified finite bounds");
    for (int axis = 0; axis < 3; ++axis)
      if (bounds_[0][axis] > (*certified)[0][axis] ||
          bounds_[1][axis] < (*certified)[1][axis])
        throw std::invalid_argument("SDF bounds do not enclose the subtree");
  }
  double Value(const Vector& point) const override { return child_->Value(point); }
  std::optional<::Box> ConservativeBounds() const override { return bounds_; }
  double LowerBound(const Vector& point) const override
  {
    if (bounds_.Inside(point))
      return -std::numeric_limits<double>::infinity();
    double squared = 0.0;
    for (int axis = 0; axis < 3; ++axis)
    {
      const double delta = point[axis] < bounds_[0][axis] ? bounds_[0][axis] - point[axis]
        : point[axis] > bounds_[1][axis] ? point[axis] - bounds_[1][axis] : 0.0;
      squared += delta * delta;
    }
    return std::sqrt(squared);
  }
private:
  NodePtr child_;
  ::Box bounds_;
};

// Common ownership and input validation for two-child construction nodes.
class BinaryNode : public Node
{
public:
  BinaryNode(NodePtr left, NodePtr right)
    : left_(std::move(left)), right_(std::move(right))
  {
    if (!left_ || !right_)
      throw std::invalid_argument("A binary SDF node requires two children");
  }

protected:
  const Node& Left() const { return *left_; }
  const Node& Right() const { return *right_; }

private:
  NodePtr left_;
  NodePtr right_;
};

} // namespace sdf
