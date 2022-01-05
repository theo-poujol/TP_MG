#include "hole_filling.h"
#include <QString>
#include <QDebug>



using namespace std ;
using namespace MeshReconstruction;

/* ************** Implicit RBF ************** */

Implicit_RBF::Implicit_RBF(vector<float> alpha, vector<float> beta, vector<MyMesh::Point> center) : _alpha(alpha), _beta(beta), _center(center)
{
    _n = _alpha.size() ;
    _d = _beta.size() ;
    if (_center.size() != _n)
        throw runtime_error("Inconsistent size of alpha and centers in Implicit_RBF constructor.");
}

// Computation of the value of the implicit surface at point X
float Implicit_RBF::val(MyMesh::Point X) const
{
    float res = 0 ;
    // Computation of the sum of RBF at centers
    for(int i=0; i<_n; i++)
    {
        res += _alpha.at(i) * myphi((X-_center.at(i)).norm()) ;
    }
    // Computation of the polynomial part
    for(int j=0; j<_d; j++)
    {
        res += _beta.at(j) * myp(j, X) ;
    }
    return res ;
}

/* ************** Hole Filling ************** */

Hole_Filling::Hole_Filling(MyMesh &mesh) : _mesh(mesh)
{
    _mesh.add_property(_vprop,  "vprop_flag");
}

// ***** Computation of boundary and its neighborhood

MyMesh::HalfedgeHandle Hole_Filling::find_boundary_edge()
{
    MyMesh::HalfedgeIter he_it = _mesh.halfedges_begin();
    while ( (he_it != _mesh.halfedges_end()) && (!_mesh.is_boundary(*he_it)))
    {
        ++he_it ;
    }
    if (he_it != _mesh.halfedges_end())
        return *he_it ;
    else
        throw std::runtime_error("Boundary HE does not exist") ;
}

std::vector<MyMesh::HalfedgeHandle> Hole_Filling::find_boundary_edges() {
    std::vector<MyMesh::HalfedgeHandle> holes;
    std::set<MyMesh::HalfedgeHandle> edge_to_ignore;

    MyMesh::HalfedgeIter he_it = _mesh.halfedges_begin();

    while ( he_it != _mesh.halfedges_end() )
    {
        if (_mesh.is_boundary(*he_it) && edge_to_ignore.find(he_it) == edge_to_ignore.end()) {
            holes.push_back(*he_it);

            MyMesh::HalfedgeHandle heh_ini = he_it ;
            MyMesh::HalfedgeHandle heh = he_it;

            // Follow (and memorize) boundary edges
            do
            {
                heh = _mesh.next_halfedge_handle(heh);
                edge_to_ignore.insert(heh);
            } while (heh != heh_ini) ;
        }

        ++he_it ;
    }

    std::cout << "Il y a " << holes.size() << " trou(s)" << std::endl;
    return holes;
}

std::vector<Hole> Hole_Filling::find_holes() {

    std::vector<Hole> holes;

    std::set<MyMesh::HalfedgeHandle> edge_to_ignore;

    MyMesh::HalfedgeIter he_it = _mesh.halfedges_begin();

    while ( he_it != _mesh.halfedges_end() )
    {
        if (_mesh.is_boundary(*he_it) && edge_to_ignore.find(he_it) == edge_to_ignore.end()) {

            std::vector<MyMesh::HalfedgeHandle> a_hole;
            a_hole.push_back(*he_it);

            MyMesh::HalfedgeHandle heh_ini = he_it ;
            MyMesh::HalfedgeHandle heh = he_it;

            // Follow (and memorize) boundary edges
            do
            {
                heh = _mesh.next_halfedge_handle(heh);
                edge_to_ignore.insert(heh);
                a_hole.push_back(heh);
            } while (heh != heh_ini) ;

            Hole hole = Hole(a_hole, _mesh);
            holes.push_back(hole);
        }

        ++he_it ;
    }

    std::vector<Hole> good_holes;
    for (Hole hole : holes) {
        if (hole.is_vertices_valid()) {
            good_holes.push_back(hole);
        }
    }

    return good_holes;
}

vector<MyMesh::VertexHandle> Hole_Filling::find_boundary(MyMesh::HalfedgeHandle heh)
{
    MyMesh::HalfedgeHandle heh_ini = heh ;
    vector<MyMesh::VertexHandle> boundary ;

    // Follow (and memorize) boundary edges
    do
    {
        boundary.push_back(_mesh.to_vertex_handle(heh));
        heh = _mesh.next_halfedge_handle(heh);
    } while (heh != heh_ini) ;
    return boundary ;
}

void Hole_Filling::init_mark_boundary(const vector<MyMesh::VertexHandle> & bnd)
{
    for (MyMesh::VertexIter v_it = _mesh.vertices_begin() ; v_it != _mesh.vertices_end(); ++v_it)
        _mesh.property(_vprop, *v_it) = false ;
    for (int i=0; i<bnd.size(); ++i)
        _mesh.property(_vprop, bnd.at(i)) = true ;
}

vector<MyMesh::VertexHandle> Hole_Filling::next_neighbors(const vector<MyMesh::VertexHandle> & bnd)
{
    // Visit bnd vertices to find and mark next circle
    vector<MyMesh::VertexHandle> next_bnd ;
    for (int i=0; i<bnd.size(); i++)
    {
        for (MyMesh::VertexVertexIter vv_it = _mesh.vv_iter(bnd.at(i));vv_it.is_valid();++vv_it)
        {
            if (_mesh.property(_vprop, *vv_it) == false) // new vertex
            {
                _mesh.property(_vprop, *vv_it) = true ;
                next_bnd.push_back(*vv_it);
            }
        }
    }
    return next_bnd ;
}

// ***** Computation of RBF

pair<pair<Eigen::MatrixXd &, Eigen::VectorXd &>, vector<MyMesh::Point> &> Hole_Filling::compute_approx_mat(vector<MyMesh::VertexHandle> vlist)
{
    const int n(vlist.size()), d(10) ;

    Eigen::MatrixXd & A = *(new Eigen::MatrixXd(3*n+d,3*n+d));
    Eigen::MatrixXd Phi(3*n,3*n);
    Eigen::MatrixXd P(3*n,d);
    Eigen::VectorXd & B = *(new Eigen::VectorXd(3*n+d));

    vector<MyMesh::Point> & pts_list = *(new vector<MyMesh::Point>);

    //Append vertices to pts_list
    for (int i=0; i<n; i++)
    {
        pts_list.push_back(_mesh.point(vlist.at(i))) ;
    }

    //Append vertices+normals to pts_list
    for (int i=0; i<n; i++)
    {
        pts_list.push_back(_mesh.point(vlist.at(i)) + _mesh.normal(vlist.at(i))) ;
    }

    //Append vertices-normals to pts_list
    for (int i=0; i<n; i++)
    {
        pts_list.push_back(_mesh.point(vlist.at(i)) - _mesh.normal(vlist.at(i))) ;
    }

    int nn = pts_list.size() ;
    // Compute corresponding B vector (0 / 1 / -1 )
    B << Eigen::VectorXd::Zero(n), Eigen::VectorXd::Ones(n), -Eigen::VectorXd::Ones(n), Eigen::VectorXd::Zero(d) ;

    // Fill Phi matrix
    // TODO

    for ( int j = 0; j < 3*n; j++) {

//      RÃ©cupÃ©ration du centre
        MyMesh::Point c = pts_list.at(j);

        for ( int i = 0; i < 3*n; i++) {

//          RÃ©cupÃ©ration du point courant
            MyMesh::Point xi = pts_list.at(i);

//          Calcul de la distance encludiÃ¨nne entre les deux points
//            double de = euclidian_distance(p2, p);
//            double e = 3.0;

//            std::cout << "de : " << de << std::endl;

            // Utilisation d'une GaussiÃ¨nne
//            Phi(j, i) = exp( - pow(( e * de), 2) );

            Phi(j,i) = myphi((xi - c).norm());
        }
    }

    // Fill P matrix
    //TODO

    for (int j = 0; j < 3*n; j++) {
        for (int i = 0; i < d; i++) {
            P(j, i) = myp(i, pts_list[j]);
        }
    }

    // Set final A matrix
    /* A = Phi | P
     *      P' | 0    */

    // TODO

    A << Phi, P, P.transpose(), Eigen::MatrixXd::Zero(d, d);

    cout << "size of pts_list : " << nn << endl ;
    cout << "End computation of matrices" << endl ;

    return {{A,B},pts_list} ;

}

pair<vector<float>&, vector<float>&> Hole_Filling::solve_approx(const pair<Eigen::MatrixXd &, Eigen::VectorXd &> &p, int n, int d)
{
    Eigen::MatrixXd & A = p.first ;
    Eigen::VectorXd & B = p.second ;

    Eigen::VectorXd res = A.householderQr().solve(B) ;

    cout << "End of solver" << endl ;
    cout << "res : " << res.head(10) << endl ;
    vector<float> & alpha = *(new vector<float>);
    vector<float> & beta = *(new vector<float>);

    if (res.size() != (n+d))
    {
        cout << "taille du res : " << res.size() << endl ;
        throw std::runtime_error("Error in solve_approx") ;
    }
    for (int i=0; i<n; i++)
    {
        alpha.push_back(res(i)) ;
    }

    for (int j=0; j<d; j++)
    {
        beta.push_back(res(n+j)) ;
    }

    return {alpha, beta} ;
}

// ***** Hole filling

// Remplir un trou
void Hole_Filling::fill_hole(string out)
{
    // TODO !!!
    MyMesh::HalfedgeHandle heh = find_boundary_edge();
    vector<MyMesh::VertexHandle> vec_vh = find_boundary(heh);

    init_mark_boundary(vec_vh);

    vector<MyMesh::VertexHandle> vec_vh_neighbors = next_neighbors(vec_vh);
    auto [approx, pts_list] = compute_approx_mat(vec_vh_neighbors);
    auto [alpha, beta] = solve_approx(approx, pts_list.size(), 10);

    Implicit_RBF rbf(alpha, beta, pts_list);
    colorize_prop();
    colorize_verts(vec_vh_neighbors);
    auto bb = estimate_BB(vec_vh_neighbors) ;

    poly_n_out(rbf, bb, out);
}

// Remplir un trou, en fonction d'un demi-arrête passée en paramètre
// idx : numéro du trou/patch, utile pour générer un fichier numéroté des patchs
void Hole_Filling::fill_a_hole(MyMesh::HalfedgeHandle heh, int idx, QString out) {
    vector<MyMesh::VertexHandle> vec_vh = find_boundary(heh);

    init_mark_boundary(vec_vh);
    vector<MyMesh::VertexHandle> vec_vh_neighbors = next_neighbors(vec_vh);
    auto [approx, pts_list] = compute_approx_mat(vec_vh_neighbors);
    auto [alpha, beta] = solve_approx(approx, pts_list.size(), 10);

    Implicit_RBF rbf(alpha, beta, pts_list);
    colorize_prop();
    colorize_verts(vec_vh_neighbors);
    auto bb = estimate_BB(vec_vh_neighbors) ;

    std::string new_out = to_out_with_idx_file_name(idx, out);

    poly_n_out(rbf, bb, new_out);
}

void Hole_Filling::fill_holes(QString out) {
    std::vector<HalfedgeHandle> holes = find_boundary_edges();

    int i = 1;
    for (auto heh : holes) {
        fill_a_hole(heh, i, out);
        i++;
    }
}


// ***** IO

void Hole_Filling::colorize_prop()
{
    for (MyMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end() ; ++v_it)
    {
        if(_mesh.property(_vprop, *v_it) == true)
            _mesh.set_color(*v_it, MyMesh::Color(255, 0, 0)) ;
        else
           _mesh.set_color(*v_it, MyMesh::Color(200, 200, 200)) ;
    }
}

void Hole_Filling::colorize_verts(const vector<MyMesh::VertexHandle> &vlist)
{
    for (MyMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end() ; ++v_it)
    {
           _mesh.set_color(*v_it, MyMesh::Color(200, 200, 200)) ;
    }

    for(int i=0; i<vlist.size(); i++)
    {
        _mesh.set_color(vlist.at(i), MyMesh::Color(255, 0, 0)) ;
    }
}

Rect3 Hole_Filling::estimate_BB(const vector<MyMesh::VertexHandle> &vlist)
{
    MyMesh::Point minp = _mesh.point(vlist.at(0)), maxp=minp ;
    for (int i=0; i<vlist.size(); i++)
    {
        minp = min(minp, _mesh.point(vlist.at(i))) ;
        maxp = max(maxp, _mesh.point(vlist.at(i))) ;
    }
    MyMesh::Point sizep(maxp-minp), centerp = (minp+maxp)/2 ;
    minp = centerp - _scale*sizep/2 ;
    sizep *= _scale ;

    Rect3 domain ;
    domain.min = {minp[0], minp[1], minp[2]} ;
    domain.size = {sizep[0], sizep[1], sizep[2]} ;
    return domain ;
}

Mesh Hole_Filling::poly_n_out(const Implicit_RBF &implicit, Rect3 domain, string filename)
{
    auto implicitEq = [&implicit](Vec3 const& pos)
    {
        MyMesh::Point X(pos.x, pos.y, pos.z) ;
        return implicit.val(X) ;
    };


    Vec3 cubeSize(domain.size*(_discr)) ;
    cout << "mind " << domain.min << endl ;
    cout << "sized" << domain.size << endl ;
    cout << "cubesize "<< cubeSize << endl ;

    auto mesh = MarchCube(implicitEq, domain, cubeSize);
    WriteObjFile(mesh, filename);

    return mesh;
}

std::string Hole_Filling::to_out_file_name(QString filename) {
    std::string out = filename.toUtf8().constData();
    unsigned int ext_pos = out.find_last_of(".");
    out = out.substr(0, ext_pos);
    out = out + "_patch.obj";

    std::cout << out << std::endl;
    return out;
}


// Génère un fichier obj du nom fileName + le numéro du patch + _patch.obj
std::string Hole_Filling::to_out_with_idx_file_name(int idx, QString filename) {
    std::string out = filename.toUtf8().constData();
    unsigned int ext_pos = out.find_last_of(".");
    out = out.substr(0, ext_pos);
    out = out + std::to_string(idx) + "_patch.obj";

    std::cout << out << std::endl;
    return out;
}
