#include "qte.h"
#include "implicits.h"
#include "sdf_field.h"
#include "sdf_operations.h"
#include "sdf_primitives.h"
#include "ui_interface.h"

#include <QtWidgets/qfiledialog.h>
#include <QtWidgets/qmessagebox.h>
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QCoreApplication>

#include <memory>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <map>

namespace {
Camera SdfViewCamera(const Vector& eye)
{
  const Vector up = std::abs(eye[2]) > 0.95 * Norm(eye)
    ? Vector::Y : Vector::Z;
  return Camera(eye, Vector(0.0), up);
}
}

MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
    uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
    connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BoxMeshExample()));
    connect(uiw->sphereImplicit, SIGNAL(clicked()), this, SLOT(SphereImplicitExample()));
    connect(uiw->generateSdf, SIGNAL(clicked()), this, SLOT(RenderSdfExample()));
    connect(uiw->exportSdfPng, SIGNAL(clicked()), this, SLOT(ExportSdfImage()));
    auto updateErosionControls = [this](int shape) {
      const bool enabled = shape == 9 || shape == 10;
      uiw->erosionCount->setEnabled(enabled);
      uiw->erosionRadius->setEnabled(enabled);
      uiw->erosionCountLabel->setEnabled(enabled);
      uiw->erosionRadiusLabel->setEnabled(enabled);
      uiw->erosionDepth->setEnabled(enabled);
      uiw->erosionSeed->setEnabled(enabled);
      uiw->erosionDepthLabel->setEnabled(enabled);
      uiw->erosionSeedLabel->setEnabled(enabled);
    };
    connect(uiw->sdfShape, QOverload<int>::of(&QComboBox::currentIndexChanged),
            this, updateErosionControls);
    updateErosionControls(uiw->sdfShape->currentIndex());
    connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
    connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
    connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BoxMeshExample()
{
	lastSdfRoot_.reset();
	Mesh boxMesh = Mesh(Box(1.0));

	std::vector<Color> cols;
	cols.resize(boxMesh.Vertexes());
    for (size_t i = 0; i < cols.size(); i++)
		cols[i] = Color(double(i) / 6.0, fmod(double(i) * 39.478378, 1.0), 0.0);

	meshColor = MeshColor(boxMesh, cols, boxMesh.VertexIndexes());
	UpdateGeometry();
}

void MainWindow::SphereImplicitExample()
{
  lastSdfRoot_.reset();
  AnalyticScalarField implicit;

  Mesh implicitMesh;
  implicit.Polygonize(31, implicitMesh, Box(2.0));

  std::vector<Color> cols;
  cols.resize(implicitMesh.Vertexes());
  for (size_t i = 0; i < cols.size(); i++)
    cols[i] = Color(0.8, 0.8, 0.8);

  meshColor = MeshColor(implicitMesh, cols, implicitMesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::RenderSdfExample()
{
  const auto started = std::chrono::steady_clock::now();

  sdf::NodePtr root;
  Mesh mesh;
  bool incrementalRendered = false;
  switch (uiw->sdfShape->currentIndex())
  {
  case 0:
    root = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
    break;
  case 1:
    root = std::make_shared<sdf::Capsule>(Vector(0.0, -0.8, 0.0),
                                          Vector(0.0, 0.8, 0.0), 0.45);
    break;
  case 2:
    root = std::make_shared<sdf::Torus>(Vector(0.0), 1.0, 0.35);
    break;
  case 3:
    root = std::make_shared<sdf::Box>(Vector(0.0), Vector(0.8, 0.6, 0.45));
    break;
  case 4:
    root = std::make_shared<sdf::Union>(
      std::make_shared<sdf::Box>(Vector(0.0), Vector(1.5, 1.5, 0.25)),
      std::make_shared<sdf::Sphere>(Vector(0.0), 1.0));
    break;
  case 5:
    root = std::make_shared<sdf::Intersection>(
      std::make_shared<sdf::Sphere>(Vector(-0.5, 0.0, 0.0), 1.1),
      std::make_shared<sdf::Sphere>(Vector(0.5, 0.0, 0.0), 1.1));
    break;
  case 6:
    root = std::make_shared<sdf::Blend>(
      std::make_shared<sdf::Sphere>(Vector(-0.8, 0.0, 0.0), 0.9),
      std::make_shared<sdf::Sphere>(Vector(0.8, 0.0, 0.0), 0.9), 0.7);
    break;
  case 7:
    root = std::make_shared<sdf::Union>(
      std::make_shared<sdf::Sphere>(Vector(-0.7, 0.0, 0.0), 0.9),
      std::make_shared<sdf::Sphere>(Vector(0.7, 0.0, 0.0), 0.9));
    break;
  case 8:
    root = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
    if (!sdf::ApplyProjectedSphericalErosion(
          root, Ray(Vector(0.0, 0.0, 3.0), Vector(0.0, 0.0, -1.0)), 0.4))
      return;
    break;
  case 9:
  case 10:
  {
    const sdf::NodePtr original = std::make_shared<sdf::Sphere>(Vector(0.0), 1.0);
    const std::vector<Ray> rays = sdf::CreateErosionProjectionRays(
      uiw->erosionCount->value(), erosionProjection_, uiw->erosionSeed->value());
    std::vector<sdf::SphericalErosionStamp> stamps =
      sdf::ProjectSphericalErosionStamps(*original, rays, uiw->erosionRadius->value());
    sdf::VarySphericalErosionRadii(stamps, uiw->erosionDepth->value(),
                                   uiw->erosionSeed->value());
    if (stamps.size() != rays.size())
      return;

    if (uiw->sdfShape->currentIndex() == 9)
    {
      uiw->generateSdf->setEnabled(false);
      meshWidget->SetCamera(SdfViewCamera(sdfCamera_));
      sdf::NodePtr removedRegion;
      for (size_t i = 0; i < stamps.size(); ++i)
      {
        removedRegion = sdf::AddSphericalErosionStamp(removedRegion, stamps[i]);
        root = sdf::MakeSphericalErosionFromUnion(original, removedRegion);
        root->Polygonize(uiw->sdfResolution->value(), mesh, Box(2.5));
        std::vector<Color> colors(mesh.Vertexes(), Color(0.8, 0.8, 0.8));
        meshColor = MeshColor(mesh, colors, mesh.VertexIndexes());
        UpdateGeometry();
        uiw->sdfGenerationTime->setText(
          QString("Impact %1 / %2").arg(i + 1).arg(stamps.size()));
        QCoreApplication::processEvents();
      }
      uiw->generateSdf->setEnabled(true);
      incrementalRendered = true;
    }
    else
    {
      root = sdf::MakeBatchedSphericalErosion(original, stamps);
    }
    break;
  }
  default:
    return;
  }

  if (!incrementalRendered)
  {
    root->Polygonize(uiw->sdfResolution->value(), mesh, Box(2.5));
    std::vector<Color> colors(mesh.Vertexes(), Color(0.8, 0.8, 0.8));
    meshColor = MeshColor(mesh, colors, mesh.VertexIndexes());
    UpdateGeometry();
  }
  meshWidget->SetCamera(SdfViewCamera(sdfCamera_));
  lastSdfRoot_ = root;
  const double milliseconds = std::chrono::duration<double, std::milli>(
    std::chrono::steady_clock::now() - started).count();
  uiw->sdfGenerationTime->setText(QString("%1 ms").arg(milliseconds, 0, 'f', 1));
}

void MainWindow::ConfigureSdfExample(int shape, int resolution)
{
  uiw->sdfShape->setCurrentIndex(shape);
  uiw->sdfResolution->setValue(resolution);
  RenderSdfExample();
}

void MainWindow::ConfigureSdfCapture(int shape, int resolution, bool wireframe,
                                     int impacts, double radius, double depth,
                                     int seed, const Vector& camera,
                                     const Vector& projection)
{
  sdfCamera_ = camera;
  erosionProjection_ = projection;
  uiw->radioShadingButton_2->setChecked(true);
  uiw->wireframe->setChecked(wireframe);
  uiw->erosionCount->setValue(impacts);
  uiw->erosionRadius->setValue(radius);
  uiw->erosionDepth->setValue(depth);
  uiw->erosionSeed->setValue(seed);
  ConfigureSdfExample(shape, resolution);
}

bool MainWindow::SaveSdfImage(const QString& path)
{
  return meshWidget->grabFramebuffer().save(path, "PNG");
}

bool MainWindow::AppendSdfStatistics(const QString& path)
{
  if (!lastSdfRoot_)
    return false;
  double maxResidual = 0.0;
  // Some clipping vertices are retained in the position buffer but never
  // referenced by a triangle; only rendered surface vertices are measured.
  std::vector<bool> referenced(meshColor.Vertexes(), false);
  for (int t = 0; t < meshColor.Triangles(); ++t)
    for (int e = 0; e < 3; ++e)
      referenced[meshColor.VertexIndex(t, e)] = true;
  for (int i = 0; i < meshColor.Vertexes(); ++i)
    if (referenced[i])
      maxResidual = std::max(maxResidual, std::abs(lastSdfRoot_->Value(meshColor.Vertex(i))));
  std::map<std::pair<int, int>, int> edgeCounts;
  for (int t = 0; t < meshColor.Triangles(); ++t)
    for (int e = 0; e < 3; ++e)
    {
      int a = meshColor.VertexIndex(t, e);
      int b = meshColor.VertexIndex(t, (e + 1) % 3);
      if (a > b) std::swap(a, b);
      ++edgeCounts[{a, b}];
    }
  int boundaryEdges = 0, nonmanifoldEdges = 0;
  for (const auto& edge : edgeCounts)
  {
    boundaryEdges += edge.second == 1;
    nonmanifoldEdges += edge.second > 2;
  }
  QFile file(path);
  const bool needsHeader = !file.exists() || file.size() == 0;
  if (!file.open(QIODevice::WriteOnly | QIODevice::Append | QIODevice::Text))
    return false;

  QTextStream stream(&file);
  if (needsHeader)
    stream << "shape_index;resolution;vertices;triangles;max_abs_sdf;boundary_edges;nonmanifold_edges\n";
  stream << uiw->sdfShape->currentIndex() << ';'
         << uiw->sdfResolution->value() << ';'
         << uiw->lineEdit->text() << ';'
         << uiw->lineEdit_2->text() << ';'
         << QString::number(maxResidual, 'g', 9) << ';'
         << boundaryEdges << ';' << nonmanifoldEdges << '\n';
  return stream.status() == QTextStream::Ok;
}

void MainWindow::ExportSdfImage()
{
  const QString suggested = QString("sdf_%1_n%2.png")
    .arg(uiw->sdfShape->currentText().toLower())
    .arg(uiw->sdfResolution->value());
  const QString path = QFileDialog::getSaveFileName(
    this, tr("Exporter la vue"), suggested, tr("Image PNG (*.png)"));
  if (path.isEmpty())
    return;

  if (!SaveSdfImage(path))
    QMessageBox::warning(this, tr("Échec de l’export"), tr("Impossible d’enregistrer l’image PNG."));
}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(lastSdfRoot_
    ? SdfViewCamera(sdfCamera_)
    : Camera(Vector(10.0, 0.0, 0.0), Vector(0.0)));
}
