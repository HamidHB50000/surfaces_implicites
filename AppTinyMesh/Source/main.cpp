#include "qte.h"
#include <QtWidgets/qapplication.h>
#include <QtCore/qtimer.h>
#include "sdf_operations.h"
#include "sdf_primitives.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

bool RunSdfBenchmark(int iterations, int repetitions, const char* outputPath)
{
  if (iterations <= 0 || repetitions <= 0 || outputPath == nullptr)
    return false;

  std::vector<Vector> points;
  points.reserve(4096);
  std::uint32_t state = 0x13579BDFu;
  for (int i = 0; i < 4096; ++i)
  {
    const auto nextCoordinate = [&state]() {
      state = 1664525u * state + 1013904223u;
      return (double(state >> 8) / double(0x00FFFFFFu)) * 4.0 - 2.0;
    };
    points.emplace_back(nextCoordinate(), nextCoordinate(), nextCoordinate());
  }

  const auto sphereA = std::make_shared<sdf::Sphere>(Vector(-0.5, 0.0, 0.0), 1.1);
  const auto sphereB = std::make_shared<sdf::Sphere>(Vector(0.5, 0.0, 0.0), 1.1);
  const auto sphereC = std::make_shared<sdf::Sphere>(Vector(0.0, 0.4, 0.0), 0.8);
  const auto box = std::make_shared<sdf::Box>(Vector(0.0), Vector(1.2, 0.9, 0.7));
  const auto torus = std::make_shared<sdf::Torus>(Vector(0.0), 1.0, 0.3);

  const std::vector<std::pair<std::string, sdf::NodePtr>> fields = {
    {"sphere", std::make_shared<sdf::Sphere>(Vector(0.0), 1.0)},
    {"capsule", std::make_shared<sdf::Capsule>(Vector(0.0, -0.8, 0.0), Vector(0.0, 0.8, 0.0), 0.45)},
    {"torus", std::make_shared<sdf::Torus>(Vector(0.0), 1.0, 0.35)},
    {"box", std::make_shared<sdf::Box>(Vector(0.0), Vector(0.8, 0.6, 0.45))},
    {"union_spheres", std::make_shared<sdf::Union>(sphereA, sphereB)},
    {"intersection_spheres", std::make_shared<sdf::Intersection>(sphereA, sphereB)},
    {"blend_spheres", std::make_shared<sdf::Blend>(sphereA, sphereB, 0.7)},
    {"nested_union_intersection", std::make_shared<sdf::Union>(
      std::make_shared<sdf::Intersection>(sphereA, box), torus)}
  };

  std::ofstream output(outputPath, std::ios::trunc);
  if (!output)
    return false;
  output << "field;iterations;repetitions;median_ms;evaluations_per_second;checksum\n";
  output << std::setprecision(12);

  for (const auto& entry : fields)
  {
    const sdf::Node& field = *entry.second;
    double checksum = 0.0;
    for (int i = 0; i < 10000; ++i)
      checksum += field.Value(points[size_t(i) % points.size()]);

    std::vector<double> elapsedMs;
    elapsedMs.reserve(size_t(repetitions));
    for (int repetition = 0; repetition < repetitions; ++repetition)
    {
      double sum = 0.0;
      const auto begin = std::chrono::steady_clock::now();
      for (int i = 0; i < iterations; ++i)
        sum += field.Value(points[size_t(i) % points.size()]);
      const auto end = std::chrono::steady_clock::now();
      checksum = sum;
      elapsedMs.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
    }

    std::sort(elapsedMs.begin(), elapsedMs.end());
    const double medianMs = elapsedMs[elapsedMs.size() / 2];
    const double evaluationsPerSecond = double(iterations) * 1000.0 / medianMs;
    output << entry.first << ';' << iterations << ';' << repetitions << ';'
           << medianMs << ';' << evaluationsPerSecond << ';' << checksum << '\n';
  }
  return bool(output);
}

bool RunBoundedNodeCheck()
{
  const auto left = std::make_shared<sdf::Sphere>(Vector(-0.5, 0.0, 0.0), 1.1);
  const auto right = std::make_shared<sdf::Sphere>(Vector(0.5, 0.0, 0.0), 1.1);
  const auto bounded = std::make_shared<sdf::BoundedNode>(
    right, Box(Vector(-0.61, -1.11, -1.11), Vector(1.61, 1.11, 1.11)));
  const sdf::Union reference(left, right);
  const sdf::Union accelerated(left, bounded);
  int compared = 0;
  for (int x = -40; x <= 40; ++x)
    for (int y = -20; y <= 20; ++y)
      for (int z = -20; z <= 20; ++z)
      {
        const Vector point(0.05 * x, 0.05 * y, 0.05 * z);
        if (std::abs(reference.Value(point) - accelerated.Value(point)) > 1e-12)
          return false;
        ++compared;
      }
  const bool distantBound = bounded->LowerBound(Vector(3.0, 0.0, 0.0)) > 0.0;
  const bool interiorBound = std::isinf(bounded->LowerBound(Vector(1.0, 0.0, 0.0))) &&
                             bounded->LowerBound(Vector(1.0, 0.0, 0.0)) < 0.0;
  bool rejectsComplement = false, rejectsUndersizedBox = false;
  try
  {
    auto complement = std::make_shared<sdf::Complement>(right);
    sdf::BoundedNode invalid(complement, Box(Vector(-0.61, -1.11, -1.11),
                                            Vector(1.61, 1.11, 1.11)));
  }
  catch (const std::invalid_argument&) { rejectsComplement = true; }
  try
  {
    sdf::BoundedNode invalid(right, Box(Vector(0.0, -1.0, -1.0),
                                       Vector(1.0, 1.0, 1.0)));
  }
  catch (const std::invalid_argument&) { rejectsUndersizedBox = true; }
  const std::vector<Ray> rays = {
    Ray(Vector(0.0, 0.0, 3.0), Vector(0.0, 0.0, -1.0)),
    Ray(Vector(-3.0, 0.2, 0.1), Vector(1.0, 0.0, 0.0)),
    Ray(Vector(0.0, 3.0, 0.0), Vector(0.0, -1.0, 0.0)),
    Ray(Vector(3.0, 3.0, 3.0), Vector(1.0, 0.0, 0.0))};
  bool sameRays = true;
  for (const Ray& ray : rays)
  {
    double plainT = 0.0, acceleratedT = 0.0;
    const bool plainHit = reference.Intersect(ray, plainT);
    const bool acceleratedHit = accelerated.Intersect(ray, acceleratedT);
    sameRays = sameRays && plainHit == acceleratedHit &&
      (!plainHit || std::abs(plainT - acceleratedT) < 1e-8);
  }
  std::cout << "Compared values: " << compared << "; bounds: "
            << distantBound << ", " << interiorBound << "; rays: " << sameRays
            << "; rejected unsafe bounds: " << rejectsComplement << ", "
            << rejectsUndersizedBox << '\n';
  return distantBound && interiorBound && sameRays && rejectsComplement &&
         rejectsUndersizedBox;
}

bool RunBoundedQueryBenchmark(int iterations, int repetitions, const char* outputPath)
{
  if (iterations <= 0 || repetitions <= 0 || outputPath == nullptr) return false;
  auto left = std::make_shared<sdf::Sphere>(Vector(-0.5, 0.0, 0.0), 1.1);
  auto right = std::make_shared<sdf::Sphere>(Vector(0.5, 0.0, 0.0), 1.1);
  auto boundedRight = std::make_shared<sdf::BoundedNode>(right,
    Box(Vector(-0.61, -1.11, -1.11), Vector(1.61, 1.11, 1.11)));
  const sdf::Union plain(left, right), accelerated(left, boundedRight);
  sdf::NodePtr heavyLeft, heavyRight;
  for (int i = 0; i < 32; ++i)
  {
    const double y = 0.28 * double(i % 4) - 0.42;
    const double z = 0.28 * double((i / 4) % 4) - 0.42;
    const double x = 0.15 * double(i / 16) - 0.25;
    const auto a = std::make_shared<sdf::Sphere>(Vector(x, y, z), 0.35);
    const auto b = std::make_shared<sdf::Sphere>(Vector(x + 6.0, y, z), 0.35);
    heavyLeft = heavyLeft ? sdf::NodePtr(std::make_shared<sdf::Union>(heavyLeft, a))
                          : sdf::NodePtr(a);
    heavyRight = heavyRight ? sdf::NodePtr(std::make_shared<sdf::Union>(heavyRight, b))
                            : sdf::NodePtr(b);
  }
  const sdf::Union heavyPlain(heavyLeft, heavyRight);
  const auto rightBounds = heavyRight->ConservativeBounds();
  if (!rightBounds) return false;
  const auto boundedHeavyRight = std::make_shared<sdf::BoundedNode>(heavyRight, *rightBounds);
  const sdf::Union heavyAccelerated(heavyLeft, boundedHeavyRight);
  std::vector<Vector> points;
  std::vector<Ray> rays;
  std::vector<Vector> localPoints;
  std::vector<Ray> localRays;
  std::uint32_t state = 0xA341316Cu;
  for (int i = 0; i < 4096; ++i)
  {
    auto coord = [&]() { state = 1664525u * state + 1013904223u;
      return 8.0 * double(state >> 8) / double(0x00FFFFFFu) - 4.0; };
    points.emplace_back(coord(), coord(), coord());
    Vector origin(coord(), coord(), coord());
    Vector target(coord() * 0.25, coord() * 0.25, coord() * 0.25);
    rays.emplace_back(origin, Normalized(target - origin));
    const Vector localPoint(0.75 * coord(), 0.4 * coord(), 0.4 * coord());
    localPoints.push_back(localPoint);
    const Vector localOrigin(0.75 * coord(), 0.5 * coord(), 0.5 * coord());
    localRays.emplace_back(localOrigin, Normalized(Vector(0.0) - localOrigin));
  }
  std::ofstream out(outputPath, std::ios::trunc);
  if (!out) return false;
  out << "scene;query;mode;iterations;repetitions;median_ms;checksum\n" << std::setprecision(12);
  for (int scene = 0; scene < 2; ++scene)
  for (int query = 0; query < 2; ++query)
    for (int mode = 0; mode < 2; ++mode)
    {
      const sdf::Node& field = scene
        ? (mode ? static_cast<const sdf::Node&>(heavyAccelerated)
                : static_cast<const sdf::Node&>(heavyPlain))
        : (mode ? static_cast<const sdf::Node&>(accelerated)
                : static_cast<const sdf::Node&>(plain));
      std::vector<double> times;
      double checksum = 0.0;
      for (int rep = 0; rep < repetitions; ++rep)
      {
        double sum = 0.0;
        auto begin = std::chrono::steady_clock::now();
        for (int i = 0; i < iterations; ++i)
        {
          const size_t index = size_t(i) % points.size();
          if (query == 0) sum += field.Value(scene ? localPoints[index] : points[index]);
          else { double t = 0.0; sum += field.Intersect(scene ? localRays[index] : rays[index], t) ? t : -1.0; }
        }
        auto end = std::chrono::steady_clock::now();
        checksum = sum;
        times.push_back(std::chrono::duration<double, std::milli>(end-begin).count());
      }
      std::sort(times.begin(), times.end());
      out << (scene ? "64_spheres" : "2_spheres") << ';'
          << (query ? "ray_intersection" : "distance") << ';'
          << (mode ? "bounded" : "plain") << ';' << iterations << ';' << repetitions
          << ';' << times[times.size()/2] << ';' << checksum << '\n';
    }
  return bool(out);
}

bool RunErosionBenchmark(int resolution, int repetitions, const char* outputPath)
{
  if (resolution < 16 || resolution > 160 || repetitions <= 0 || outputPath == nullptr)
    return false;

  const sdf::NodePtr original = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
  const std::vector<Ray> rays = sdf::CreateErosionProjectionRays(64, Vector(5.0, 4.0, 6.0));
  std::vector<sdf::SphericalErosionStamp> allStamps =
    sdf::ProjectSphericalErosionStamps(*original, rays, 0.12);
  sdf::VarySphericalErosionRadii(allStamps);
  if (allStamps.size() != rays.size())
    return false;

  std::ofstream output(outputPath, std::ios::trunc);
  if (!output)
    return false;
  output << "stamp_count;mode;resolution;repetitions;median_total_ms;vertices;triangles\n";
  output << std::setprecision(12);

  for (int stampCount : {8, 32, 64})
  {
    const std::vector<sdf::SphericalErosionStamp> stamps(
      allStamps.begin(), allStamps.begin() + stampCount);
    for (int mode = 0; mode < 2; ++mode)
    {
      std::vector<double> elapsedMs;
      elapsedMs.reserve(size_t(repetitions));
      int vertices = 0;
      int triangles = 0;
      for (int repetition = 0; repetition < repetitions; ++repetition)
      {
        const auto begin = std::chrono::steady_clock::now();
        if (mode == 0)
        {
          sdf::NodePtr removedRegion;
          sdf::NodePtr incremental;
          Mesh mesh;
          for (const sdf::SphericalErosionStamp& stamp : stamps)
          {
            removedRegion = sdf::AddSphericalErosionStamp(removedRegion, stamp);
            incremental = sdf::MakeSphericalErosionFromUnion(original, removedRegion);
            incremental->Polygonize(resolution, mesh, Box(2.5));
          }
          vertices = mesh.Vertexes();
          triangles = mesh.Triangles();
        }
        else
        {
          sdf::NodePtr batch = sdf::MakeBatchedSphericalErosion(original, stamps);
          Mesh mesh;
          batch->Polygonize(resolution, mesh, Box(2.5));
          vertices = mesh.Vertexes();
          triangles = mesh.Triangles();
        }
        const auto end = std::chrono::steady_clock::now();
        elapsedMs.push_back(std::chrono::duration<double, std::milli>(end - begin).count());
      }

      std::sort(elapsedMs.begin(), elapsedMs.end());
      const double medianMs = elapsedMs[elapsedMs.size() / 2];
      output << stampCount << ';' << (mode == 0 ? "incremental" : "batch") << ';'
             << resolution << ';' << repetitions << ';' << medianMs << ';'
             << vertices << ';' << triangles << '\n';
    }
  }
  return bool(output);
}

} // namespace

int main(int argc, char *argv[])
{
	if (argc == 2 && std::string(argv[1]) == "--test-bounded-node")
		return RunBoundedNodeCheck() ? 0 : 1;
	if (argc == 5 && std::string(argv[1]) == "--benchmark-bounded")
	{
		try
		{
			return RunBoundedQueryBenchmark(std::stoi(argv[2]), std::stoi(argv[3]), argv[4]) ? 0 : 4;
		}
		catch (const std::exception& error)
		{
			std::cerr << "Invalid bounded benchmark arguments: " << error.what() << '\n';
			return 2;
		}
	}
	if (argc == 5 && std::string(argv[1]) == "--benchmark-sdf")
	{
		try
		{
			const int iterations = std::stoi(argv[2]);
			const int repetitions = std::stoi(argv[3]);
			return RunSdfBenchmark(iterations, repetitions, argv[4]) ? 0 : 4;
		}
		catch (const std::exception& error)
		{
			std::cerr << "Invalid benchmark arguments: " << error.what() << '\n';
			return 2;
		}
	}
	if (argc == 5 && std::string(argv[1]) == "--benchmark-erosion")
	{
		try
		{
			const int resolution = std::stoi(argv[2]);
			const int repetitions = std::stoi(argv[3]);
			return RunErosionBenchmark(resolution, repetitions, argv[4]) ? 0 : 4;
		}
		catch (const std::exception& error)
		{
			std::cerr << "Invalid erosion benchmark arguments: " << error.what() << '\n';
			return 2;
		}
	}

	QApplication app(argc, argv);

	MainWindow mainWin;
	mainWin.showMaximized();

	// Save the rendered OpenGL view and optional mesh statistics.
	if (argc >= 5 && std::string(argv[1]) == "--capture-sdf")
	{
		try
		{
			const int shape = std::stoi(argv[2]);
			const int resolution = std::stoi(argv[3]);
			const QString path = QString::fromLocal8Bit(argv[4]);
			QString statisticsPath;
			bool wireframe = false;
			int impacts = 64, seed = 0;
			double radius = 0.12, depth = 0.55;
			Vector camera(5.0, 4.0, 6.0);
			Vector projection(5.0, 4.0, 6.0);
			for (int i = 5; i < argc; ++i)
			{
				const std::string option(argv[i]);
				if (option == "--wireframe") wireframe = true;
				else if (option == "--impacts" && i + 1 < argc) impacts = std::stoi(argv[++i]);
				else if (option == "--radius" && i + 1 < argc) radius = std::stod(argv[++i]);
				else if (option == "--depth" && i + 1 < argc) depth = std::stod(argv[++i]);
				else if (option == "--seed" && i + 1 < argc) seed = std::stoi(argv[++i]);
				else if (option == "--camera" && i + 3 < argc)
				{
					const double x = std::stod(argv[++i]);
					const double y = std::stod(argv[++i]);
					const double z = std::stod(argv[++i]);
					camera = Vector(x, y, z);
				}
				else if (option == "--projection" && i + 3 < argc)
				{
					const double x = std::stod(argv[++i]);
					const double y = std::stod(argv[++i]);
					const double z = std::stod(argv[++i]);
					projection = Vector(x, y, z);
				}
				else if (option.rfind("--", 0) != 0 && statisticsPath.isEmpty())
					statisticsPath = QString::fromLocal8Bit(argv[i]);
				else throw std::invalid_argument("Unknown or incomplete capture option");
			}
			if (shape < 0 || shape > 10 || resolution < 16 || resolution > 160 ||
				impacts < 1 || impacts > 96 || radius < 0.02 || radius > 0.40 ||
				depth < 0.0 || depth > 1.0 || seed < 0 || seed > 9999 ||
				!std::isfinite(radius) || !std::isfinite(depth) ||
				!std::isfinite(Norm(camera)) || Norm(camera) < 1e-12 ||
				!std::isfinite(Norm(projection)) || Norm(projection) < 1e-12)
				return 2;

			QTimer::singleShot(500, &mainWin, [&mainWin, shape, resolution, path,
				statisticsPath, wireframe, impacts, radius, depth, seed, camera, projection]() {
				mainWin.ConfigureSdfCapture(shape, resolution, wireframe, impacts,
					radius, depth, seed, camera, projection);
				QTimer::singleShot(500, &mainWin, [&mainWin, path, statisticsPath]() {
					const bool savedImage = mainWin.SaveSdfImage(path);
					const bool savedStatistics = statisticsPath.isEmpty() ||
						mainWin.AppendSdfStatistics(statisticsPath);
					QApplication::exit(savedImage && savedStatistics ? 0 : 3);
				});
			});
		}
		catch (const std::exception& error)
		{
			std::cerr << "Invalid capture arguments: " << error.what() << '\n';
			return 2;
		}
	}

	return app.exec();
}
