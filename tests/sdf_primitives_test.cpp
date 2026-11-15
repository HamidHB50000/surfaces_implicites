#include "sdf_primitives.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>

namespace
{

void Near(double actual, double expected)
{
  assert(std::abs(actual - expected) < 1e-12);
}

template <typename F>
void Rejects(F makeInvalid)
{
  bool rejected = false;
  try { makeInvalid(); }
  catch (const std::invalid_argument&) { rejected = true; }
  assert(rejected);
}

} // namespace

int main()
{
  const sdf::Sphere sphere(Vector(1.0, 2.0, 3.0), 2.0);
  Near(sphere.Value(Vector(1.0, 2.0, 3.0)), -2.0);
  Near(sphere.Value(Vector(3.0, 2.0, 3.0)), 0.0);
  Near(sphere.Value(Vector(4.0, 2.0, 3.0)), 1.0);

  const sdf::Capsule capsule(Vector(0.0, 0.0, 0.0), Vector(0.0, 2.0, 0.0), 0.5);
  Near(capsule.Value(Vector(0.0, 1.0, 0.0)), -0.5);
  Near(capsule.Value(Vector(0.5, 1.0, 0.0)), 0.0);
  Near(capsule.Value(Vector(0.0, 3.0, 0.0)), 0.5);
  const sdf::Capsule collapsed(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 0.0), 0.5);
  Near(collapsed.Value(Vector(1.0, 0.0, 0.0)), 0.5);

  const sdf::Torus torus(Vector(0.0, 0.0, 0.0), 2.0, 0.5);
  Near(torus.Value(Vector(2.0, 0.0, 0.0)), -0.5);
  Near(torus.Value(Vector(2.5, 0.0, 0.0)), 0.0);
  Near(torus.Value(Vector(0.0, 0.0, 0.0)), 1.5);

  const sdf::Box box(Vector(1.0, 0.0, 0.0), Vector(1.0, 2.0, 3.0));
  Near(box.Value(Vector(1.0, 0.0, 0.0)), -1.0);
  Near(box.Value(Vector(2.0, 0.0, 0.0)), 0.0);
  Near(box.Value(Vector(3.0, 3.0, 0.0)), std::sqrt(2.0));
  Near(box.Value(Vector(2.0, 2.0, 3.0)), 0.0);

  const sdf::NodePtr polymorphic = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
  Near(polymorphic->Value(Vector(2.0, 0.0, 0.0)), 1.0);

  Rejects([] { sdf::Sphere(Vector(0.0), 0.0); });
  Rejects([] { sdf::Torus(Vector(0.0), 0.5, 0.5); });
  Rejects([] { sdf::Box(Vector(0.0), Vector(1.0, -1.0, 1.0)); });
  Rejects([] { sdf::Sphere(Vector(std::numeric_limits<double>::infinity(), 0.0, 0.0), 1.0); });
}
