#include "hole.h"

#include "edge.h"

class BowyerWatson {
public :

    BowyerWatson(std::vector<Point> points) : points(points), super_triangle(SuperTriangle(Point(0,0),Point(0,0),Point(0,0))) {
        std::vector<Point> super_points = sort_points();
        assert(super_points.size() == 3);
        super_triangle = SuperTriangle(super_points.at(0), super_points.at(1), super_points.at(2));
    }

    SuperTriangle super_triangle;

    std::vector<Point> points;

    std::vector<Point> sort_points() {
        std::vector<Point> super_points;

        float xmin = 0, xmax = 0, ymin = 0, ymax = 0;


        for (Point p : this->points) {
            if (p.x < xmin) xmin = p.x;
            if (p.y < ymin) ymin = p.y;

            if (p.x > xmax) xmax = p.x;
            if (p.y > ymax) ymax = p.y;
        }

        Point * A = new Point ( ((xmax+xmin)/2), ymax);
        Point * B = new Point ( xmin, ymin);
        Point * C = new Point ( xmax, ymin);

//        Point * A = new Point ( xmin - 2, ymin - 2);
//        Point * B = new Point ( xmin + 2*(xmax - xmin) + 3*2, ymin - 2);
//        Point * C = new Point ( xmin - 2, ymin + 2 * (ymax - ymin) + 3*2);

        super_points.push_back(*A);
        super_points.push_back(*B);
        super_points.push_back(*C);

        return super_points;
    }

    // Triangulation de Delaunay
    std::vector<Triangle> Delaunay() {

        std::set<Triangle> fait;
        std::vector<Triangle> triangles;

        triangles.push_back(this->super_triangle);

        int i = 0;
        for (Point point : this->points) {

            std::vector<Triangle> englobant;
            std::vector<Triangle> a_supprimer;

            for (Triangle triangle : triangles) {
                Point * centre = triangle.circum;
                float rayon = triangle.circum_radius;

                if (triangle.isInside(point)) {
                    i++;
                    englobant.push_back(triangle);
                    a_supprimer.push_back(triangle);
                }
            }

            for (Triangle triangle : a_supprimer) {
                auto it = std::find(triangles.begin(), triangles.end(), triangle);
                if (it != triangles.end()) triangles.erase(it);
            }

            for (Triangle triangle : englobant) {
                std::set<Triangle> test;
                Triangle * AB_Point = new Triangle(triangle.A, triangle.B, point);
                Triangle * BC_Point = new Triangle(triangle.B, triangle.C, point);
                Triangle * CA_Point = new Triangle(triangle.C, triangle.A, point);

                triangles.push_back(*AB_Point);
                triangles.push_back(*BC_Point);
                triangles.push_back(*CA_Point);
            }
        }

    std::cout << "i : " << i << std::endl;
    return triangles;
    }
};

Hole::Hole(vector<MyMesh::HalfedgeHandle> hl, MyMesh &mesh) : _mesh(mesh)
{
    this->he_list = hl;

    std::set<EdgeHandle> ignore_edges;

    std::set<VertexHandle> ignore_vertices;

    for (auto heh : this->he_list) {
        EdgeHandle eh = this->_mesh.edge_handle(heh);

        if (ignore_edges.find(eh) == ignore_edges.end()) {
            this->e_list.push_back(eh);
            ignore_edges.insert(eh);
        }
    }

    for (auto eh : this->e_list) {

        auto vh1 = mesh.to_vertex_handle(mesh.halfedge_handle(eh, 0));
        auto vh2 = mesh.from_vertex_handle(mesh.halfedge_handle(eh, 0));

        if (ignore_vertices.find(vh1) == ignore_vertices.end()) {
            this->pts_list.push_back(vh1);
            ignore_vertices.insert(vh1);
        }

        if (ignore_vertices.find(vh2) == ignore_vertices.end()) {
            this->pts_list.push_back(vh2);
            ignore_vertices.insert(vh2);
        }
    }


}

bool Hole::is_valid() {

    if (this->e_list.size() < 3) return false;

    return true;
}

bool Hole::is_vertices_valid() {

    for (auto vh : this->pts_list) {
        int cpt_face = 0;
        for (MyMesh::VertexFaceIter vf_iter = _mesh.vf_iter(vh); vf_iter.is_valid(); vf_iter++) {
            cpt_face++;
        }

        // Si le sommet en question n'a qu'une seule face adjacente
        if (cpt_face == 1) return false;
    }

    return true;
}

void Hole::Bowyer_Watson() {
    assert(this->pts_list.size() >= 3);

    std::vector<Point> points;

    for (int i = 0; i < this->pts_list.size(); i++) {
        MyMesh::Point myPoint = _mesh.point(this->pts_list.at(i));
        Point * p = new Point(myPoint[0], myPoint[1]);
        p->setVh(this->pts_list.at(i));
        points.push_back(*p);  
    }

    BowyerWatson * bw = new BowyerWatson(points);
    std::vector<Triangle> triangles = bw->Delaunay();

    std::vector<Triangle> bons;

    std::cout << "trou a " << this->pts_list.size() << " sommets" << std::endl;
    std::cout << "triangles : " << triangles.size() << std::endl;

    for (Triangle t : triangles) {
        auto itA = std::find(points.begin(), points.end(), t.A);
        auto itB = std::find(points.begin(), points.end(), t.B);
        auto itC = std::find(points.begin(), points.end(), t.C);

        if (itA != points.end() && itB != points.end() && itC != points.end()) bons.push_back(t);

    }

    std::cout << "il y a " << bons.size() << " bons triangles" << std::endl;

    for (Triangle t : bons) {

        _mesh.set_color(*t.A.vh, MyMesh::Color(0, 0, 255 ));
        _mesh.data(*t.A.vh).thickness = 9;

        _mesh.set_color(*t.B.vh, MyMesh::Color(0, 0, 255 ));
        _mesh.data(*t.B.vh).thickness = 9;

        _mesh.set_color(*t.C.vh, MyMesh::Color(0, 0, 255 ));
        _mesh.data(*t.C.vh).thickness = 9;

        _mesh.add_face(*t.A.vh, *t.B.vh, *t.C.vh);
    }
}

void Hole::fix() {
    int nb_vertices = this->pts_list.size();

    std::vector<std::vector<VertexHandle>> vertices;
    if (nb_vertices == 3) {
        _mesh.add_face(this->pts_list);
    }

    else {
        this->Bowyer_Watson();
    }
}
