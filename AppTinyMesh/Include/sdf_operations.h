#pragma once

#include "sdf_node.h"

#include <vector>

namespace sdf
{

// Negates a field while reversing its inside/outside convention.
class Complement final : public Node
{
public:
  explicit Complement(NodePtr child);
  double Value(const Vector& point) const override;

private:
  NodePtr child_;
};

class Union final : public BinaryNode
{
public:
  using BinaryNode::BinaryNode;
  double Value(const Vector& point) const override;
  void Polygonize(int resolution, Mesh& mesh, const ::Box& bounds,
                  const double& epsilon = 1e-4) const override;
  std::optional<::Box> ConservativeBounds() const override;
};

class Intersection final : public BinaryNode
{
public:
  using BinaryNode::BinaryNode;
  double Value(const Vector& point) const override;
  void Polygonize(int resolution, Mesh& mesh, const ::Box& bounds,
                  const double& epsilon = 1e-4) const override;
};

// Polynomial smooth union. Width controls the distance range of the blend.
class Blend final : public BinaryNode
{
public:
  Blend(NodePtr left, NodePtr right, double width);
  double Value(const Vector& point) const override;

private:
  double width_;
};

struct SphericalErosionStamp
{
  Vector center;
  double radius;
  Vector outwardNormal;
};

// Projects a spherical stamp onto a field without modifying it.
bool ProjectSphericalErosionStamp(const AnalyticScalarField& field, const Ray& ray,
                                  double radius, SphericalErosionStamp& stamp,
                                  double* hitParameter = nullptr);
std::vector<Ray> CreateErosionProjectionRays(int count, const Vector& viewingDirection,
                                             int seed = 0);
std::vector<SphericalErosionStamp> ProjectSphericalErosionStamps(
  const AnalyticScalarField& field, const std::vector<Ray>& rays, double radius);
void VarySphericalErosionRadii(std::vector<SphericalErosionStamp>& stamps,
                               double depth = 0.55, int seed = 0);

// Applies one already-projected stamp to an object field.
void ApplySphericalErosion(NodePtr& object, const SphericalErosionStamp& stamp);

// Incrementally appends one stamp to the removed-region union.
NodePtr AddSphericalErosionStamp(NodePtr removedRegion,
                                 const SphericalErosionStamp& stamp);
NodePtr MakeSphericalErosionFromUnion(NodePtr object, NodePtr removedRegion);

// Applies a set of stamps in one balanced union and one subtraction.
NodePtr MakeBatchedSphericalErosion(NodePtr object,
                                    const std::vector<SphericalErosionStamp>& stamps);

// Projects a spherical brush along a ray and subtracts it from the object.
// Returns false and leaves object unchanged when the ray misses.
bool ApplyProjectedSphericalErosion(NodePtr& object, const Ray& projectionRay,
                                   double radius, double* hitParameter = nullptr);

} // namespace sdf
