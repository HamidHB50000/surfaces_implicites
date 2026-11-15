#pragma once

#include "sdf_node.h"

namespace sdf
{

class Sphere final : public Node
{
public:
  Sphere(const Vector& center, double radius);
  double Value(const Vector& point) const override;
  std::optional<::Box> ConservativeBounds() const override;

private:
  Vector center_;
  double radius_;
};

class Capsule final : public Node
{
public:
  Capsule(const Vector& start, const Vector& end, double radius);
  double Value(const Vector& point) const override;
  std::optional<::Box> ConservativeBounds() const override;

private:
  Vector start_;
  Vector end_;
  double radius_;
};

// Ring torus centered at center, with its symmetry axis parallel to Y.
class Torus final : public Node
{
public:
  Torus(const Vector& center, double majorRadius, double minorRadius);
  double Value(const Vector& point) const override;
  std::optional<::Box> ConservativeBounds() const override;

private:
  Vector center_;
  double majorRadius_;
  double minorRadius_;
};

// Axis-aligned box. Half-extents are measured from its center.
class Box final : public Node
{
public:
  Box(const Vector& center, const Vector& halfExtents);
  double Value(const Vector& point) const override;
  std::optional<::Box> ConservativeBounds() const override;

private:
  Vector center_;
  Vector halfExtents_;
};

} // namespace sdf
