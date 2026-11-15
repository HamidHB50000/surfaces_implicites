#include "sdf_operations.h"
#include "sdf_primitives.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

int main()
{
  sdf::NodePtr object = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
  double hitParameter = -1.0;
  const bool eroded = sdf::ApplyProjectedSphericalErosion(
    object, Ray(Vector(0.0, 0.0, 3.0), Vector(0.0, 0.0, -1.0)), 0.4, &hitParameter);

  assert(eroded);
  assert(std::abs(hitParameter - 2.0) < 2e-4);
  assert(object->Value(Vector(0.0, 0.0, 1.0)) > 0.0); // Removed by the brush.
  assert(object->Value(Vector(0.0, 0.0, 0.4)) < 0.0); // Material remains below the crater.
  assert(object->Value(Vector(0.8, 0.0, 0.5)) < 0.0); // Distant surface remains intact.

  const sdf::Node* unchanged = object.get();
  const bool missed = sdf::ApplyProjectedSphericalErosion(
    object, Ray(Vector(3.0, 0.0, 3.0), Vector(0.0, 0.0, -1.0)), 0.4);
  assert(!missed);
  assert(object.get() == unchanged);

  bool rejectedRadius = false;
  try
  {
    sdf::ApplyProjectedSphericalErosion(
      object, Ray(Vector(0.0, 0.0, 3.0), Vector(0.0, 0.0, -1.0)), 0.0);
  }
  catch (const std::invalid_argument&)
  {
    rejectedRadius = true;
  }
  assert(rejectedRadius);

  // Incremental updates and one balanced batch must define the same field.
  const sdf::NodePtr original = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
  const std::vector<Ray> rays = sdf::CreateErosionProjectionRays(64, Vector(5.0, 4.0, 6.0));
  std::vector<sdf::SphericalErosionStamp> stamps =
    sdf::ProjectSphericalErosionStamps(*original, rays, 0.12);
  sdf::VarySphericalErosionRadii(stamps);
  assert(stamps.size() == rays.size());

  sdf::NodePtr removedRegion;
  sdf::NodePtr incremental;
  for (const sdf::SphericalErosionStamp& stamp : stamps)
  {
    removedRegion = sdf::AddSphericalErosionStamp(removedRegion, stamp);
    incremental = sdf::MakeSphericalErosionFromUnion(original, removedRegion);
  }
  const sdf::NodePtr batch = sdf::MakeBatchedSphericalErosion(original, stamps);

  for (int x = -6; x <= 6; ++x)
    for (int y = -6; y <= 6; ++y)
      for (int z = -6; z <= 6; ++z)
      {
        const Vector point(0.25 * x, 0.25 * y, 0.25 * z);
        assert(std::abs(incremental->Value(point) - batch->Value(point)) < 1e-12);
      }
}
