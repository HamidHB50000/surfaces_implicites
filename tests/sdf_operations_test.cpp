#include "sdf_operations.h"
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
  const sdf::NodePtr left = std::make_shared<sdf::Sphere>(Vector(-1.0, 0.0, 0.0), 1.0);
  const sdf::NodePtr right = std::make_shared<sdf::Sphere>(Vector(1.0, 0.0, 0.0), 1.0);
  const sdf::Union united(left, right);
  const sdf::Intersection common(left, right);
  const sdf::Blend blended(left, right, 1.0);

  Near(united.Value(Vector(0.0)), 0.0);
  Near(common.Value(Vector(0.0)), 0.0);
  Near(blended.Value(Vector(0.0)), -0.25);
  Near(united.Value(Vector(-1.0, 0.0, 0.0)), -1.0);
  Near(common.Value(Vector(-1.0, 0.0, 0.0)), 1.0);
  Near(blended.Value(Vector(3.0, 0.0, 0.0)), 1.0);

  // A composed tree retains shared leaves and evaluates recursively.
  const sdf::NodePtr joined = std::make_shared<sdf::Union>(left, right);
  const sdf::NodePtr box = std::make_shared<sdf::Box>(Vector(0.0), Vector(0.5));
  const sdf::Intersection tree(joined, box);
  Near(tree.Value(Vector(0.0)), 0.0);
  Near(tree.Value(Vector(0.5, 0.0, 0.0)), 0.0);
  Near(tree.Value(Vector(3.0, 0.0, 0.0)), 2.5);

  Rejects([&] { sdf::Union bad(nullptr, right); });
  Rejects([&] { sdf::Intersection bad(left, nullptr); });
  Rejects([&] { sdf::Blend bad(left, right, 0.0); });
  Rejects([&] { sdf::Blend bad(left, right, std::numeric_limits<double>::infinity()); });
}
