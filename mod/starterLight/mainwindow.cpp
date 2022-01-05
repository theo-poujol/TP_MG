#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "HoleFilling/hole_filling.h"
#include "hole.h"


/* **** début de la partie boutons et IHM **** */

vector<vector<VertexHandle>> MainWindow::detectCrack(MyMesh *_mesh){
    vector<vector<VertexHandle>> boundaries;
    EdgeHandle eh;
    int i = 0;


    // parcours de toutes les aretes
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        // trouvé une arrete qui a une seule face adjacente
        HalfedgeHandle heh0 = _mesh->halfedge_handle(*curEdge,0);
        if(!_mesh->face_handle(heh0).is_valid() || !_mesh->opposite_face_handle(heh0).is_valid() ){
            HalfedgeHandle iter;
            HalfedgeHandle init;
            //récupère l'halfedges en bordure
            if(!_mesh->face_handle(heh0).is_valid()){
                iter = heh0;
                init = heh0;
            }
            if(!_mesh->opposite_face_handle(heh0).is_valid()){
                iter = _mesh->opposite_halfedge_handle(heh0);
                init = _mesh->opposite_halfedge_handle(heh0);
            }
            iter = _mesh->next_halfedge_handle(iter);
            if(_mesh->data(_mesh->to_vertex_handle(iter)).taken == true)
                continue;
            boundaries.push_back(vector<VertexHandle>());
            boundaries[i].push_back(_mesh->to_vertex_handle(iter));
            _mesh->data(_mesh->to_vertex_handle(iter)).taken = true;
            //edge length
            holes_length.push_back(_mesh->calc_edge_length(*curEdge));
            //iter sur la bordure
            while(iter != init) {
                iter = _mesh->next_halfedge_handle(iter);
                //cout << _mesh->next_halfedge_handle(iter) << endl;
                holes_length[i] += _mesh->calc_edge_length(_mesh->edge_handle(iter));
                boundaries[i].push_back(_mesh->to_vertex_handle(iter));
                _mesh->data(_mesh->to_vertex_handle(iter)).taken = true;
            }
            i++;
            //break;
        }
    }
    return boundaries;
}


std::vector<FaceHandle> MainWindow::identifyOverlapes() {

    bool isDeleted[mesh.n_faces()];
    std::vector<FaceHandle> overlaped_faces;

    for (MyMesh::FaceIter curFace = mesh.faces_begin(); curFace != mesh.faces_end(); curFace++) {
        for (MyMesh::FaceVertexIter curVer = mesh.fv_iter(*curFace); curVer.is_valid(); curVer++) {
            VertexHandle vh = *curVer;

            int cpt_face = 0;

            for (MyMesh::VertexFaceIter vf_iter = mesh.vf_iter(vh); vf_iter.is_valid(); vf_iter++) {
                cpt_face++;
            }

            if (cpt_face == 1) {
                FaceHandle fh = *curFace;
                if (!isDeleted[fh.idx()]) {
                    overlaped_faces.push_back(fh);
                    isDeleted[fh.idx()] = true;
                }

            }
        }
    }

    return overlaped_faces;
}

int MainWindow::find_closest_vertex(VertexHandle vh)
{
    float min = 10000;
    Vec3f pointVh = mesh.point(vh);
    //VertexHandle closestVer;
    int closestVer;
    for(int i = 0; i< (int)cracks[1].size(); i++){
        Vec3f pointCur = mesh.point(cracks[1][i]);
        if(norm(pointVh - pointCur) < min){
            min = norm(pointVh - pointCur);
            closestVer = i;
            continue;
        }
    }
    return closestVer;
}


// exemple pour charger un fichier .obj
void MainWindow::on_pushButton_chargement_clicked()
{
    // fenêtre de sélection des fichiers
    QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "", tr("Mesh Files (*.obj)"));

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
}

// Identification des trous
void MainWindow::on_pushButton_identify_holes_clicked()
{

    resetAllColorsAndThickness(&mesh);

    Hole_Filling hf(mesh);

//    std::vector<MyMesh::HalfedgeHandle> holes = hf.find_boundary_edges();

    std::vector<Hole> holes = hf.find_holes();

    for (Hole hole : holes) {

        if (hole.is_vertices_valid()) {
            for (auto eh : hole.e_list) {

                mesh.set_color(eh, MyMesh::Color(0, 255, 0));
                mesh.data(eh).thickness = 7;
            }

            for (auto vh : hole.pts_list) {
                mesh.set_color(vh, MyMesh::Color(255, 0,0 ));
                mesh.data(vh).thickness = 9;
            }
        }
    }

    displayMesh(&mesh);
}

void MainWindow::on_pushButton_fix_holes_clicked() {
    resetAllColorsAndThickness(&mesh);

    Hole_Filling hf(mesh);

    std::vector<Hole> holes = hf.find_holes();

    if (holes.size() == 0) std::cout << "Il n'y a aucuns trous de détecter." << std::endl;

    while (holes.size() != 0) {

        for (Hole hole : holes) {
            hole.fix();
        }
        holes = hf.find_holes();
    }

    displayMesh(&mesh);
}

void MainWindow::on_pushButton_identify_overlaping_clicked() {

    std::vector<FaceHandle> overlaped_faces = identifyOverlapes();

    if (overlaped_faces.size() > 0) {
        for (auto fh : overlaped_faces) {
            mesh.set_color(fh, MyMesh::Color(0,255,0 ));
        }
    }

    std::cout << "taille :" << overlaped_faces.size() << std::endl;

    displayMesh(&mesh);
}

void MainWindow::on_pushButton_detectFis_clicked()
{
    vector<vector<VertexHandle>> holes = detectCrack(&mesh);
    int max = 0;

    for(int lenght : holes_length){
        if(lenght > max)
            max = lenght;
    }
    for(int i = 0; i < (int)holes.size(); i++){
        int color = 255;
        if((int)holes_length[i] == max){
            cracks.push_back(holes[i]);
            for(VertexHandle vh : holes[i]){
                mesh.data(vh).thickness = 6;
                mesh.set_color(vh,MyMesh::Color(color,0,0));
            }
        }
    }
    displayMesh(&mesh);
}


/* **** fin de la partie boutons et IHM **** */


/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 1;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }

    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(150, 150, 150));
    }

    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(0, 0, 0));
    }
}

// charge un objet MyMesh dans l'environnement OpenGL
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange)
{
    GLuint* triIndiceArray = new GLuint[_mesh->n_faces() * 3];
    GLfloat* triCols = new GLfloat[_mesh->n_faces() * 3 * 3];
    GLfloat* triVerts = new GLfloat[_mesh->n_faces() * 3 * 3];

    int i = 0;

    if(isTemperatureMap)
    {
        QVector<float> values;

        if(mapRange == -1)
        {
            for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                values.append(fabs(_mesh->data(*curVert).value));
            qSort(values);
            mapRange = values.at(values.size()*0.8);
            qDebug() << "mapRange" << mapRange;
        }

        float range = mapRange;
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;

        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            if(_mesh->data(*fvIt).value > 0){triCols[3*i+0] = 255; triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            else{triCols[3*i+2] = 255; triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0); triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);}
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }
    else
    {
        MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
        MyMesh::ConstFaceVertexIter fvIt;
        for (; fIt!=fEnd; ++fIt)
        {
            fvIt = _mesh->cfv_iter(*fIt);
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++; ++fvIt;
            triCols[3*i+0] = _mesh->color(*fIt)[0]; triCols[3*i+1] = _mesh->color(*fIt)[1]; triCols[3*i+2] = _mesh->color(*fIt)[2];
            triVerts[3*i+0] = _mesh->point(*fvIt)[0]; triVerts[3*i+1] = _mesh->point(*fvIt)[1]; triVerts[3*i+2] = _mesh->point(*fvIt)[2];
            triIndiceArray[i] = i;

            i++;
        }
    }


    ui->displayWidget->loadMesh(triVerts, triCols, _mesh->n_faces() * 3 * 3, triIndiceArray, _mesh->n_faces() * 3);

    delete[] triIndiceArray;
    delete[] triCols;
    delete[] triVerts;

    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit)
    {
        float t = _mesh->data(*eit).thickness;
        if(t > 0)
        {
            if(!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext())
    {
        it.next();

        for(int e = 0; e < it.value().size(); e++)
        {
            int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh1)[0];
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;

    GLuint* pointsIndiceArray = new GLuint[_mesh->n_vertices()];
    GLfloat* pointsCols = new GLfloat[_mesh->n_vertices() * 3];
    GLfloat* pointsVerts = new GLfloat[_mesh->n_vertices() * 3];

    i = 0;
    QHash<float, QList<int> > vertsIDbyThickness;
    for (MyMesh::VertexIter vit = _mesh->vertices_begin(); vit != _mesh->vertices_end(); ++vit)
    {
        float t = _mesh->data(*vit).thickness;
        if(t > 0)
        {
            if(!vertsIDbyThickness.contains(t))
                vertsIDbyThickness[t] = QList<int>();
            vertsIDbyThickness[t].append((*vit).idx());
        }
    }
    QHashIterator<float, QList<int> > vitt(vertsIDbyThickness);
    QList<QPair<float, int> > vertsSizes;

    while (vitt.hasNext())
    {
        vitt.next();

        for(int v = 0; v < vitt.value().size(); v++)
        {
            int vidx = vitt.value().at(v);

            pointsVerts[3*i+0] = _mesh->point(_mesh->vertex_handle(vidx))[0];
            pointsVerts[3*i+1] = _mesh->point(_mesh->vertex_handle(vidx))[1];
            pointsVerts[3*i+2] = _mesh->point(_mesh->vertex_handle(vidx))[2];
            pointsCols[3*i+0] = _mesh->color(_mesh->vertex_handle(vidx))[0];
            pointsCols[3*i+1] = _mesh->color(_mesh->vertex_handle(vidx))[1];
            pointsCols[3*i+2] = _mesh->color(_mesh->vertex_handle(vidx))[2];
            pointsIndiceArray[i] = i;
            i++;
        }
        vertsSizes.append(qMakePair(vitt.key(), vitt.value().size()));
    }

    ui->displayWidget->loadPoints(pointsVerts, pointsCols, i * 3, pointsIndiceArray, i, vertsSizes);

    delete[] pointsIndiceArray;
    delete[] pointsCols;
    delete[] pointsVerts;
}


MainWindow::MainWindow(QWidget *parent) : QMainWindow(parent), ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_repFis_clicked()
{
    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();

    std::vector<MyMesh::VertexHandle> uneNouvelleFace;
    for (int i = 0; i < (int)cracks[0].size(); i++) {
        int closest = find_closest_vertex(cracks[0][i]);
        int next = i+1;
        if (i == (int)cracks[0].size()-1)
            next = 0;
        int closest2 = find_closest_vertex(cracks[0][next]);
        if (closest2 == closest){
            cout << "same" << endl;
            closest2 = closest-1;
        }

        uneNouvelleFace.clear();
        uneNouvelleFace.push_back(cracks[0][i]);
        uneNouvelleFace.push_back(cracks[0][next]);
        uneNouvelleFace.push_back(cracks[1][closest2]);
        uneNouvelleFace.push_back(cracks[1][closest]);
        mesh.add_face(uneNouvelleFace);
    }

    for(int i = 0; i < (int)cracks[0].size(); i++){
        int closest = find_closest_vertex(cracks[0][i]);


        if(mesh.is_collapse_ok(mesh.find_halfedge(cracks[0][i], cracks[1][closest])))
            mesh.collapse(mesh.find_halfedge(cracks[0][i], cracks[1][closest]));

    }

    mesh.garbage_collection();

    resetAllColorsAndThickness(&mesh);

    displayMesh(&mesh);
}

void MainWindow::on_pushButton_fix_overlaping_clicked()
{

    mesh.request_vertex_status();
    mesh.request_edge_status();
    mesh.request_face_status();
    int cpt = 0;

    std::vector<FaceHandle> faces = identifyOverlapes();

    std::cout << "s : " << faces.size() << std::endl;

    if (faces.size() > 0) {
        for (auto fh : faces ) {
            mesh.delete_face(fh, true);
        }
    }

    mesh.garbage_collection();

    resetAllColorsAndThickness(&mesh);

    displayMesh(&mesh);
}
