#include "sdf_operations.h"
#include "sdf_primitives.h"

#include <cassert>
#include <cmath>
#include <memory>

namespace
{

void IntersectsAt(const sdf::Node& field, const Ray& ray, double expectedT,
                  double tTolerance = 2e-4)
{
  double t = -1.0;
  assert(field.Intersect(ray, t));
  assert(std::abs(t - expectedT) < tTolerance);
  assert(std::abs(field.Value(ray(t))) < 2e-4);
}

void Misses(const sdf::Node& field, const Ray& ray)
{
  double t = -1.0;
  assert(!field.Intersect(ray, t));
  assert(t == 0.0);
}

} // namespace

int main()
{
  const sdf::Sphere sphere(Vector(0.0), 1.0);

  // First hit from outside, including conversion back to a non-unit ray parameter.
  IntersectsAt(sphere, Ray(Vector(-3.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0)), 2.0);
  IntersectsAt(sphere, Ray(Vector(-3.0, 0.0, 0.0), Vector(2.0, 0.0, 0.0)), 1.0);

  // A tangent ray converges to the surface; a nearby ray must remain a miss.
  IntersectsAt(sphere, Ray(Vector(-2.0, 1.0, 0.0), Vector(1.0, 0.0, 0.0)), 2.0, 0.02);
  Misses(sphere, Ray(Vector(-2.0, 1.01, 0.0), Vector(1.0, 0.0, 0.0)));

  // Surface and interior starts return the current hit and the forward exit.
  IntersectsAt(sphere, Ray(Vector(1.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0)), 0.0);
  IntersectsAt(sphere, Ray(Vector(0.0), Vector(1.0, 0.0, 0.0)), 1.0);

  // A zero direction is invalid and must not report a hit.
  Misses(sphere, Ray(Vector(-3.0, 0.0, 0.0), Vector(0.0)));

  // The inherited intersection routine also evaluates a composed SDF tree.
  const sdf::NodePtr left = std::make_shared<sdf::Sphere>(Vector(-1.0, 0.0, 0.0), 1.0);
  const sdf::NodePtr right = std::make_shared<sdf::Sphere>(Vector(1.0, 0.0, 0.0), 1.0);
  const sdf::Union united(left, right);
  IntersectsAt(united, Ray(Vector(-4.0, 0.0, 0.0), Vector(1.0, 0.0, 0.0)), 2.0);
}
