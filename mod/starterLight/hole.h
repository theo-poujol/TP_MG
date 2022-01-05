#ifndef HOLE_H
#define HOLE_H
#include <iostream>

#include "mainwindow.h"
#include "supertriangle.h"

using namespace std;


class Hole
{
public:
    Hole(vector<MyMesh::HalfedgeHandle> hl, MyMesh &mesh);
    MyMesh &_mesh ;
    vector<MyMesh::HalfedgeHandle> he_list;
    vector<MyMesh::EdgeHandle> e_list;
    vector<MyMesh::VertexHandle> pts_list;
    vector<MyMesh::FaceHandle> face_list;


    bool is_valid();
    bool is_vertices_valid();

    void Bowyer_Watson();
    void fix();
};

#endif // HOLE_H
