#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QFileDialog>
#include <QMainWindow>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace Ui {
class MainWindow;
}

using namespace OpenMesh;
using namespace OpenMesh::Attributes;

struct MyTraits : public OpenMesh::DefaultTraits
{
    // use vertex normals and vertex colors
    VertexAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    // store the previous halfedge
    HalfedgeAttributes( OpenMesh::Attributes::PrevHalfedge );
    // use face normals face colors
    FaceAttributes( OpenMesh::Attributes::Normal | OpenMesh::Attributes::Color | OpenMesh::Attributes::Status);
    EdgeAttributes( OpenMesh::Attributes::Color | OpenMesh::Attributes::Status );
    // vertex thickness
    VertexTraits{float thickness; float value; bool taken;};
    // edge thickness
    EdgeTraits{float thickness;};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MyTraits> MyMesh;

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:

    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

    void displayMesh(MyMesh *_mesh, bool isTemperatureMap = false, float mapRange = -1);
    void resetAllColorsAndThickness(MyMesh* _mesh);
    std::vector<std::vector<VertexHandle>> detectCrack(MyMesh *_mesh);
    int find_closest_vertex(VertexHandle vh);
    std::vector<FaceHandle> identifyOverlapes();

    bool inArray(std::vector<int> const& visited, int const& a);
    int findNewFace(std::vector<int> const& visited, int const& nbFaces);
    std::vector<int> exploreNeighbors(MyMesh* _mesh, int face, std::vector<int> const& visited, int const& color);
    std::vector<std::vector<std::vector<int>>> detecteParts(MyMesh* _mesh);
    void showParts(MyMesh* _mesh, std::vector<std::vector<std::vector<int>>> noises);
    bool isNoice(std::vector<int> const& parts, MyMesh* _mesh);
    std::vector<int> vertexEquivalence(MyMesh* _mesh);
    bool is_equal(MyMesh* _mesh, int v1, int v2);
    void deleteBruit(MyMesh* _mesh, std::vector<std::vector<int>> noise);

    void fixall();

private slots:
    void on_pushButton_suppBruit_clicked();
    void on_pushButton_suppBruitMaillage_clicked();
    void on_pushButton_chargement_clicked();
    void on_pushButton_identify_holes_clicked();
    void on_pushButton_fix_holes_clicked();
    void on_pushButton_identify_overlaping_clicked();
    void on_pushButton_fix_overlaping_clicked();
    void on_pushButton_detectFis_clicked();
    void on_pushButton_repFis_clicked();
    void on_pushButton_voirBruit_clicked();
    void on_pushButton_fixall_clicked();

private:

    MyMesh mesh;
    std::vector<float> holes_length;
    std::vector<std::vector<VertexHandle>> cracks;
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
