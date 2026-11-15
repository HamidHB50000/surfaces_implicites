#ifndef __Qte__
#define __Qte__

#include <QtWidgets/qmainwindow.h>
#include "realtime.h"
#include "meshcolor.h"
#include "sdf_node.h"

QT_BEGIN_NAMESPACE
	namespace Ui { class Assets; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
  Q_OBJECT
private:
  Ui::Assets* uiw;           //!< Interface

  MeshWidget* meshWidget;   //!< Viewer
  MeshColor meshColor;		//!< Mesh.
  sdf::NodePtr lastSdfRoot_;
  Vector sdfCamera_ = Vector(5.0, 4.0, 6.0);
  Vector erosionProjection_ = Vector(5.0, 4.0, 6.0);

public:
  MainWindow();
  ~MainWindow();
  void CreateActions();
  void UpdateGeometry();
  void ConfigureSdfExample(int shape, int resolution);
  void ConfigureSdfCapture(int shape, int resolution, bool wireframe,
                           int impacts, double radius, double depth, int seed,
                           const Vector& camera, const Vector& projection);
  bool SaveSdfImage(const QString& path);
  bool AppendSdfStatistics(const QString& path);

public slots:
  void editingSceneLeft(const Ray&);
  void editingSceneRight(const Ray&);
  void BoxMeshExample();
  void SphereImplicitExample();
  void RenderSdfExample();
  void ExportSdfImage();
  void ResetCamera();
  void UpdateMaterial();
};

#endif
