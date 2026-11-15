#include "sdf_primitives.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace sdf
{
namespace
{

void CheckPoint(const Vector& point)
{
  for (int i = 0; i < 3; ++i)
    if (!std::isfinite(point[i]))
      throw std::invalid_argument("SDF coordinates must be finite");
}

void CheckPositive(double value)
{
  if (!std::isfinite(value) || value <= 0.0)
    throw std::invalid_argument("SDF dimensions must be finite and positive");
}

} // namespace

Sphere::Sphere(const Vector& center, double radius)
  : center_(center), radius_(radius)
{
  CheckPoint(center_);
  CheckPositive(radius_);
}

double Sphere::Value(const Vector& point) const
{
  return Norm(point - center_) - radius_;
}

std::optional<::Box> Sphere::ConservativeBounds() const
{
  return ::Box(center_ - Vector(radius_), center_ + Vector(radius_));
}

Capsule::Capsule(const Vector& start, const Vector& end, double radius)
  : start_(start), end_(end), radius_(radius)
{
  CheckPoint(start_);
  CheckPoint(end_);
  CheckPositive(radius_);
}

double Capsule::Value(const Vector& point) const
{
  const Vector axis = end_ - start_;
  const double lengthSquared = SquaredNorm(axis);
  if (lengthSquared == 0.0)
    return Norm(point - start_) - radius_;

  const double t = std::clamp(((point - start_) * axis) / lengthSquared, 0.0, 1.0);
  return Norm(point - (start_ + axis * t)) - radius_;
}

std::optional<::Box> Capsule::ConservativeBounds() const
{
  Vector low, high;
  for (int axis = 0; axis < 3; ++axis)
  {
    low[axis] = std::min(start_[axis], end_[axis]) - radius_;
    high[axis] = std::max(start_[axis], end_[axis]) + radius_;
  }
  return ::Box(low, high);
}

Torus::Torus(const Vector& center, double majorRadius, double minorRadius)
  : center_(center), majorRadius_(majorRadius), minorRadius_(minorRadius)
{
  CheckPoint(center_);
  CheckPositive(majorRadius_);
  CheckPositive(minorRadius_);
  if (majorRadius_ <= minorRadius_)
    throw std::invalid_argument("A ring torus requires major radius > minor radius");
}

double Torus::Value(const Vector& point) const
{
  const Vector p = point - center_;
  return std::hypot(std::hypot(p[0], p[2]) - majorRadius_, p[1]) - minorRadius_;
}

std::optional<::Box> Torus::ConservativeBounds() const
{
  const Vector extent(majorRadius_ + minorRadius_, minorRadius_,
                      majorRadius_ + minorRadius_);
  return ::Box(center_ - extent, center_ + extent);
}

Box::Box(const Vector& center, const Vector& halfExtents)
  : center_(center), halfExtents_(halfExtents)
{
  CheckPoint(center_);
  CheckPoint(halfExtents_);
  for (int i = 0; i < 3; ++i)
    CheckPositive(halfExtents_[i]);
}

double Box::Value(const Vector& point) const
{
  const Vector p = point - center_;
  const double qx = std::abs(p[0]) - halfExtents_[0];
  const double qy = std::abs(p[1]) - halfExtents_[1];
  const double qz = std::abs(p[2]) - halfExtents_[2];
  const double outside = std::hypot(std::max(qx, 0.0),
                                    std::max(qy, 0.0),
                                    std::max(qz, 0.0));
  const double inside = std::min(std::max(qx, std::max(qy, qz)), 0.0);
  return outside + inside;
}

std::optional<::Box> Box::ConservativeBounds() const
{
  return ::Box(center_ - halfExtents_, center_ + halfExtents_);
}

} // namespace sdf
