#include "sdf_operations.h"
#include "sdf_primitives.h"

#include <algorithm>
#include <cmath>
#include <array>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace sdf
{
namespace
{

void PolygonizeHardBoundary(const Node& field, const Node& left, const Node& right,
                            bool leftWinsWhenLess, int resolution, Mesh& mesh,
                            const ::Box& bounds, const double& epsilon)
{
  field.AnalyticScalarField::Polygonize(resolution, mesh, bounds, epsilon);

  std::vector<Vector> positions;
  std::vector<double> difference;
  positions.reserve(mesh.Vertexes());
  difference.reserve(mesh.Vertexes());
  for (int i = 0; i < mesh.Vertexes(); ++i)
  {
    positions.push_back(mesh.Vertex(i));
    difference.push_back(left.Value(positions.back()) - right.Value(positions.back()));
  }

  std::map<std::pair<int, int>, int> crossings;
  const auto crossing = [&](int a, int b) -> int
  {
    if (difference[a] == 0.0) return a;
    if (difference[b] == 0.0) return b;

    const auto edge = std::minmax(a, b);
    const auto found = crossings.find(edge);
    if (found != crossings.end()) return found->second;

    const double t = difference[a] / (difference[a] - difference[b]);
    Vector point = positions[a] + t * (positions[b] - positions[a]);
    const Vector linearPoint = point;
    const double edgeLength = Norm(positions[b] - positions[a]);

    // Seek the intersection of the child surfaces while sharing each edge
    // vertex between adjacent triangles. Keep the interpolated point if the
    // local projection cannot be validated.
    for (int iteration = 0; iteration < 12; ++iteration)
    {
      const double fa = left.Value(point);
      const double fb = right.Value(point);
      if (std::max(std::abs(fa), std::abs(fb)) < 1e-7) break;

      const Vector ga = left.Gradient(point);
      const Vector gb = right.Gradient(point);
      const double aa = ga * ga;
      const double ab = ga * gb;
      const double bb = gb * gb;
      const double determinant = aa * bb - ab * ab;
      if (determinant < 1e-10) break;

      Vector step = ((-fa * bb + fb * ab) * ga + (-fb * aa + fa * ab) * gb)
        / determinant;
      const double stepLength = Norm(step);
      const double maxStep = 0.25 * edgeLength;
      if (stepLength > maxStep)
        step *= maxStep / stepLength;
      point += step;
    }

    const bool local = Norm(point - linearPoint) <= 0.75 * edgeLength;
    const bool onBoth = std::abs(left.Value(point)) <= 5e-3 &&
                        std::abs(right.Value(point)) <= 5e-3;
    if (!local || !onBoth || std::abs(field.Value(point)) > 5e-3)
      point = linearPoint;

    const int index = static_cast<int>(positions.size());
    positions.push_back(point);
    difference.push_back(0.0);
    crossings.emplace(edge, index);
    return index;
  };

  std::vector<int> vertexIndexes;
  std::vector<int> normalIndexes;
  const std::vector<int> triangles = mesh.VertexIndexes();
  vertexIndexes.reserve(triangles.size() * 2);
  normalIndexes.reserve(triangles.size() * 2);

  for (size_t i = 0; i + 2 < triangles.size(); i += 3)
  {
    const std::array<int, 3> original = {triangles[i], triangles[i + 1], triangles[i + 2]};
    for (int side = 0; side < 2; ++side)
    {
      std::vector<int> polygon;
      polygon.reserve(4);
      for (int j = 0; j < 3; ++j)
      {
        const int a = original[j];
        const int b = original[(j + 1) % 3];
        const bool insideA = side == 0
          ? (leftWinsWhenLess ? difference[a] <= 0.0 : difference[a] >= 0.0)
          : (leftWinsWhenLess ? difference[a] > 0.0 : difference[a] < 0.0);
        const bool insideB = side == 0
          ? (leftWinsWhenLess ? difference[b] <= 0.0 : difference[b] >= 0.0)
          : (leftWinsWhenLess ? difference[b] > 0.0 : difference[b] < 0.0);
        if (insideA) polygon.push_back(a);
        if (insideA != insideB) polygon.push_back(crossing(a, b));
      }

      if (polygon.size() > 1 && polygon.front() == polygon.back())
        polygon.pop_back();
      for (size_t j = 1; j + 1 < polygon.size(); ++j)
      {
        if (polygon[0] == polygon[j] || polygon[j] == polygon[j + 1] ||
            polygon[j + 1] == polygon[0])
          continue;

        const int ids[] = {polygon[0], polygon[j], polygon[j + 1]};
        for (const int id : ids)
        {
          vertexIndexes.push_back(id);
          normalIndexes.push_back(2 * id + side);
        }
      }
    }
  }

  std::vector<Vector> normals;
  normals.reserve(2 * positions.size());
  for (const Vector& point : positions)
  {
    normals.push_back(left.Normal(point));
    normals.push_back(right.Normal(point));
  }
  mesh = Mesh(positions, normals, vertexIndexes, normalIndexes);
}

NodePtr BuildBalancedUnion(const std::vector<NodePtr>& nodes, size_t begin, size_t end)
{
  if (end - begin == 1)
    return nodes[begin];
  const size_t middle = begin + (end - begin) / 2;
  return std::make_shared<Union>(BuildBalancedUnion(nodes, begin, middle),
                                 BuildBalancedUnion(nodes, middle, end));
}

} // namespace

Complement::Complement(NodePtr child) : child_(std::move(child))
{
  if (!child_)
    throw std::invalid_argument("A complement node requires a child");
}

double Complement::Value(const Vector& point) const
{
  return -child_->Value(point);
}

double Union::Value(const Vector& point) const
{
  const double leftValue = Left().Value(point);
  if (Right().LowerBound(point) >= leftValue)
    return leftValue;
  return std::min(leftValue, Right().Value(point));
}

std::optional<::Box> Union::ConservativeBounds() const
{
  const auto left = Left().ConservativeBounds();
  const auto right = Right().ConservativeBounds();
  if (!left || !right) return std::nullopt;
  Vector low, high;
  for (int axis = 0; axis < 3; ++axis)
  {
    low[axis] = std::min((*left)[0][axis], (*right)[0][axis]);
    high[axis] = std::max((*left)[1][axis], (*right)[1][axis]);
  }
  return ::Box(low, high);
}

void Union::Polygonize(int resolution, Mesh& mesh, const ::Box& bounds,
                       const double& epsilon) const
{
  PolygonizeHardBoundary(*this, Left(), Right(), true, resolution, mesh, bounds, epsilon);
}

double Intersection::Value(const Vector& point) const
{
  return std::max(Left().Value(point), Right().Value(point));
}

void Intersection::Polygonize(int resolution, Mesh& mesh, const ::Box& bounds,
                              const double& epsilon) const
{
  PolygonizeHardBoundary(*this, Left(), Right(), false, resolution, mesh, bounds, epsilon);
}

Blend::Blend(NodePtr left, NodePtr right, double width)
  : BinaryNode(std::move(left), std::move(right)), width_(width)
{
  if (!std::isfinite(width_) || width_ <= 0.0)
    throw std::invalid_argument("SDF blend width must be finite and positive");
}

double Blend::Value(const Vector& point) const
{
  const double a = Left().Value(point);
  const double b = Right().Value(point);
  const double h = std::max(width_ - std::abs(a - b), 0.0) / width_;
  return std::min(a, b) - 0.25 * width_ * h * h;
}

bool ApplyProjectedSphericalErosion(NodePtr& object, const Ray& projectionRay,
                                    double radius, double* hitParameter)
{
  if (!object)
    throw std::invalid_argument("Spherical erosion requires an object field");
  SphericalErosionStamp stamp;
  if (!ProjectSphericalErosionStamp(*object, projectionRay, radius, stamp, hitParameter))
    return false;
  ApplySphericalErosion(object, stamp);
  return true;
}

bool ProjectSphericalErosionStamp(const AnalyticScalarField& field, const Ray& ray,
                                  double radius, SphericalErosionStamp& stamp,
                                  double* hitParameter)
{
  if (!std::isfinite(radius) || radius <= 0.0)
    throw std::invalid_argument("Erosion sphere radius must be finite and positive");

  double hit = 0.0;
  if (!field.Intersect(ray, hit))
    return false;

  const Vector center = ray(hit);
  const Vector gradient = field.Gradient(center);
  const double gradientLength = Norm(gradient);
  const Vector outward = std::isfinite(gradientLength) && gradientLength > 1e-12
    ? gradient / gradientLength
    : -ray.Direction() / Norm(ray.Direction());
  stamp = SphericalErosionStamp{center, radius, outward};
  if (hitParameter != nullptr)
    *hitParameter = hit;
  return true;
}

std::vector<Ray> CreateErosionProjectionRays(int count, const Vector& viewingDirection,
                                             int seed)
{
  if (count <= 0)
    throw std::invalid_argument("Erosion projection ray count must be positive");

  const double viewLength = Norm(viewingDirection);
  if (!std::isfinite(viewLength) || viewLength <= 1e-12)
    throw std::invalid_argument("Erosion viewing direction must be finite and nonzero");
  const Vector view = viewingDirection / viewLength;
  const Vector reference = std::abs(view[1]) < 0.9
    ? Vector(0.0, 1.0, 0.0)
    : Vector(1.0, 0.0, 0.0);
  const Vector tangent = Normalized(view / reference);
  const Vector bitangent = Normalized(tangent / view);
  std::vector<Ray> rays;
  rays.reserve(size_t(count));
  for (int i = 0; i < count; ++i)
  {
    const double x = 0.78 * std::sin(12.9898 * (i + 1) + 3.1 * std::sin(i * 2.7) + seed * 0.731);
    const double y = 0.78 * std::sin(7.233 * (i + 1) + 2.3 * std::cos(i * 1.9) + seed * 1.193);
    const Vector direction = Normalized(view + x * tangent + y * bitangent);
    rays.emplace_back(3.0 * direction, -direction);
  }
  return rays;
}

std::vector<SphericalErosionStamp> ProjectSphericalErosionStamps(
  const AnalyticScalarField& field, const std::vector<Ray>& rays, double radius)
{
  if (!std::isfinite(radius) || radius <= 0.0)
    throw std::invalid_argument("Erosion sphere radius must be finite and positive");

  std::vector<SphericalErosionStamp> stamps;
  stamps.reserve(rays.size());
  for (const Ray& ray : rays)
  {
    SphericalErosionStamp stamp;
    if (ProjectSphericalErosionStamp(field, ray, radius, stamp))
      stamps.push_back(stamp);
  }
  return stamps;
}

void VarySphericalErosionRadii(std::vector<SphericalErosionStamp>& stamps,
                               double depth, int seed)
{
  if (!std::isfinite(depth) || depth < 0.0 || depth > 1.0)
    throw std::invalid_argument("Erosion depth must be in [0, 1]");
  for (size_t i = 0; i < stamps.size(); ++i)
  {
    stamps[i].radius *= 0.55 + 0.85 * (0.5 + 0.5 * std::sin(2.9 * double(i) + 0.7 + seed * 0.417));
    const double normalLength = Norm(stamps[i].outwardNormal);
    if (!std::isfinite(normalLength) || normalLength <= 1e-12)
      throw std::invalid_argument("Erosion stamp needs a valid outward normal");
    stamps[i].center += (1.0 - depth) * stamps[i].radius *
                        (stamps[i].outwardNormal / normalLength);
  }
}

void ApplySphericalErosion(NodePtr& object, const SphericalErosionStamp& stamp)
{
  if (!object)
    throw std::invalid_argument("Spherical erosion requires an object field");
  if (!std::isfinite(stamp.radius) || stamp.radius <= 0.0)
    throw std::invalid_argument("Erosion sphere radius must be finite and positive");

  const NodePtr sphere = std::make_shared<Sphere>(stamp.center, stamp.radius);
  object = std::make_shared<Intersection>(
    object, std::make_shared<Complement>(sphere));
}

NodePtr AddSphericalErosionStamp(NodePtr removedRegion,
                                 const SphericalErosionStamp& stamp)
{
  if (!std::isfinite(stamp.radius) || stamp.radius <= 0.0)
    throw std::invalid_argument("Erosion sphere radius must be finite and positive");
  const NodePtr sphere = std::make_shared<Sphere>(stamp.center, stamp.radius);
  if (!removedRegion)
    return sphere;
  return std::make_shared<Union>(std::move(removedRegion), sphere);
}

NodePtr MakeSphericalErosionFromUnion(NodePtr object, NodePtr removedRegion)
{
  if (!object || !removedRegion)
    throw std::invalid_argument("Erosion requires an object and a nonempty removed region");
  return std::make_shared<Intersection>(
    std::move(object), std::make_shared<Complement>(std::move(removedRegion)));
}

NodePtr MakeBatchedSphericalErosion(NodePtr object,
                                    const std::vector<SphericalErosionStamp>& stamps)
{
  if (!object)
    throw std::invalid_argument("Spherical erosion requires an object field");
  if (stamps.empty())
    throw std::invalid_argument("A batched erosion requires at least one stamp");

  std::vector<NodePtr> spheres;
  spheres.reserve(stamps.size());
  for (const SphericalErosionStamp& stamp : stamps)
  {
    if (!std::isfinite(stamp.radius) || stamp.radius <= 0.0)
      throw std::invalid_argument("Erosion sphere radius must be finite and positive");
    spheres.push_back(std::make_shared<Sphere>(stamp.center, stamp.radius));
  }

  const NodePtr removedRegion = BuildBalancedUnion(spheres, 0, spheres.size());
  return MakeSphericalErosionFromUnion(std::move(object), removedRegion);
}

} // namespace sdf
