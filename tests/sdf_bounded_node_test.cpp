#include "sdf_operations.h"
#include "sdf_primitives.h"

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>

int main()
{
  const auto left = std::make_shared<sdf::Sphere>(Vector(-0.5, 0.0, 0.0), 1.1);
  const auto rawRight = std::make_shared<sdf::Sphere>(Vector(0.5, 0.0, 0.0), 1.1);
  const auto boundedRight = std::make_shared<sdf::BoundedNode>(
    rawRight, Box(Vector(-0.61, -1.11, -1.11), Vector(1.61, 1.11, 1.11)));
  const sdf::Union plain(left, rawRight);
  const sdf::Union accelerated(left, boundedRight);

  for (int x = -40; x <= 40; ++x)
    for (int y = -20; y <= 20; ++y)
      for (int z = -20; z <= 20; ++z)
      {
        const Vector p(0.05 * x, 0.05 * y, 0.05 * z);
        assert(std::abs(plain.Value(p) - accelerated.Value(p)) < 1e-12);
      }

  assert(boundedRight->LowerBound(Vector(3.0, 0.0, 0.0)) > 0.0);
  assert(std::isinf(boundedRight->LowerBound(Vector(1.0, 0.0, 0.0))));
  bool rejectedComplement = false;
  try
  {
    sdf::BoundedNode invalid(std::make_shared<sdf::Complement>(rawRight),
                             Box(Vector(-0.61, -1.11, -1.11),
                                 Vector(1.61, 1.11, 1.11)));
  }
  catch (const std::invalid_argument&) { rejectedComplement = true; }
  assert(rejectedComplement);
}
