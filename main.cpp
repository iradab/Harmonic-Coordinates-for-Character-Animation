// INF 585 Project
// Harmonic Coordinates
//  Gael Vanderlee and Irada Bunyatova

#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <iostream>
#include <igl/unproject_onto_mesh.h>
#include <ostream>
#include <chrono>
#include <thread>
#include <utility>


using namespace Eigen;
using namespace std;

class Grid {
public:
    int side_length = std::pow(2, 6);
    float cell_step_x;
    float cell_step_y;
    vector<vector<float>> cell_types;
    vector<vector<vector<float>>>cell_vectors;
    MatrixXd Vcage;
    float min_x = 0.0;
    float min_y = 0.0;

    Grid(int width, int height) {
        cell_step_x = width / side_length;
        cell_step_y = height / side_length;

        for (int j = 0; j < side_length; j++) {
            vector<float> v1;
            for (int k = 0; k < side_length; k++) {
                v1.push_back(0);
            }
            cell_types.push_back(v1);
        }
    }

    void updateTypes(MatrixXd cage_points) {
        Vcage = cage_points;
        for (int i = 0; i < side_length * side_length; i++) {
            int col = i % side_length;
            int row = i / side_length;
            cell_types[row][col] = getCellType(row, col);
        }
        for (int i = 0; i < side_length; i++) {
            vector<vector<float>> v2;
            for (int j = 0; j < side_length; j++) {
                vector<float> v1;
                for (int k = 0; k < Vcage.rows(); k++) {
                    v1.push_back(0);
                }
                v2.push_back(v1);
            }
            cell_vectors.push_back(v2);
        }
    }


    int getCellType(int row, int col) {
        RowVector3d bottom_left = RowVector3d(min_x + col * cell_step_x,
            min_y + row * cell_step_y, 0);
        RowVector3d bottom_right = RowVector3d(min_x + col * cell_step_x + cell_step_x,
            min_y + row * cell_step_y, 0);
        RowVector3d top_left = RowVector3d(min_x + col * cell_step_x,
            min_y + row * cell_step_y + cell_step_y, 0);
        RowVector3d top_right = RowVector3d(min_x + col * cell_step_x + cell_step_x,
            min_y + row * cell_step_y + cell_step_y, 0);

        bool is_in[4] = { PointInPolygon(bottom_left, Vcage),
                         PointInPolygon(bottom_right, Vcage),
                         PointInPolygon(top_left, Vcage),
                         PointInPolygon(top_right, Vcage) };


        bool inside = true;
        bool outside = true;

        for (bool i : is_in) {
            if (i)
                outside = false;
            else
                inside = false;
        }

        if (inside) return 0;
        else if (outside) return 2;
        else return 1;
    }

    bool PointInPolygon(RowVector3d point, MatrixXd polygon) {

        int i, j, nvert = polygon.rows();
        bool is_inside = false;
        for (i = 0, j = nvert - 1; i < nvert; j = i++) {

            if (((polygon(i, 1) >= point(1)) != (polygon(j, 1) >= point(1))) &&
                (point(0) <= (polygon(j, 0) - polygon(i, 0)) * (point(1) - polygon(i, 1)) /
                    (polygon(j, 1) - polygon(i, 1)) + polygon(i, 0))
                )
                is_inside = !is_inside;
        }
        return is_inside;
    }

    void updateCoordinatesAtBoundary() {
        int v1, v2;
        for (int i = 0; i < Vcage.rows(); i++) {
            if (i == Vcage.rows() - 1) {
                v1 = i;
                v2 = 0;
            }
            else {
                v1 = i;
                v2 = i + 1;
            }

            int y1 = int(Vcage.row(v1)[0] / cell_step_x);
            int x1 = int(Vcage.row(v1)[1] / cell_step_y);

            int y2 = int(Vcage.row(v2)[0] / cell_step_x);
            int x2 = int(Vcage.row(v2)[1] / cell_step_y);

            int minx = min(x1, x2), maxx = max(x1, x2);
            int miny = min(y1, y2), maxy = max(y1, y2);

            int nbBcells = 0;
            if (minx == maxx)
                nbBcells = maxy - miny + 1;
            if (miny == maxy)
                nbBcells = maxx - minx + 1;
            for (int a = minx; a <= maxx; a++) {
                for (int b = miny; b <= maxy; b++) {
                    if (cell_types[a][b] == 1)
                        nbBcells++;
                }
            }

            nbBcells = nbBcells * 10;
            for (int interpStep = 0; interpStep <= nbBcells; interpStep++) {
                RowVector3d p = interpStep * (1.f / nbBcells) * Vcage.row(v1) + (1 - interpStep * (1.f / nbBcells)) * Vcage.row(v2);
                int yp = int(p(0) / cell_step_x);
                int xp = int(p(1) / cell_step_y);

                cell_vectors[xp][yp][v1] = interpStep * (1.f / nbBcells);
                cell_vectors[xp][yp][v2] = (1 - interpStep * (1.f / nbBcells));

            }
        }
    }
    void updateCoordinatesInBoundary() {
        bool continuing = true;
        while (continuing == true) {
            continuing = false;
            for (int i = 0; i < side_length; i++) {
                for (int j = 0; j < side_length; j++) {
                    if (cell_types[i][j] == 0) {
                        float averageChange = 0;
                        for (int v = 0; v < Vcage.rows(); v++) {
                            float initial = cell_vectors[i][j][v];
                            cell_vectors[i][j][v] = (cell_vectors[i + 1][j][v] + cell_vectors[i - 1][j][v] + cell_vectors[i][j + 1][v] + cell_vectors[i][j - 1][v]) / 4;
                            averageChange += abs(initial - cell_vectors[i][j][v]);
                        }
                        if (averageChange / Vcage.rows() > 0.00001) {
                            continuing = true;
                        }
                    }
                }
            }
        }

    }
    void updateCoordinates() {
        updateCoordinatesAtBoundary();
        updateCoordinatesInBoundary();
    }

    void updateMesh(MatrixXd& Mesh, MatrixXd& Cage) {
        Vcage = Cage;
        std::cout << "Updating mesh.." << std::endl;
        for (int i = 0; i < Mesh.rows(); i++) {
            //std::cout << "Mesh.row(" << i << ") = " << Mesh.row(i) << std::endl;
            if (Mesh.row(i) != RowVector3d(0, 0, 0)) {
                float sum = 0;
                int y = int(Mesh(i, 0) / cell_step_x);
                int x = int(Mesh(i, 1) / cell_step_y);

                Mesh.row(i) = RowVector3d(0, 0, 0);

                for (int j = 0; j < Vcage.rows(); j++) {
                    Mesh.row(i) += cell_vectors[x][y][j] * Vcage.row(j);
                    sum += cell_vectors[x][y][j];
                }
                //std::cout << "Mesh2.row(" << i << ") = " << Mesh.row(i) << std::endl;
                std::cout << "Sum of harmonic coordinates for vertex n."<< i << " is = " << sum << std::endl;
            }
        }
    }
    void draw_grid(igl::opengl::glfw::Viewer& viewer) {
        std::vector<RowVector3d> Vins, Vout, Vbound;

        for (int i = 0; i < side_length; i++) {
            for (int j = 0; j < side_length; j++) {
                if (cell_types[i][j] == 0) { // for inside ones
                    Vins.push_back(RowVector3d(min_x + j * cell_step_x + cell_step_x / 2,
                        min_y + i * cell_step_y + cell_step_y / 2, 0));
                }
                if (cell_types[i][j] == 1) { // for boundary ones   // || cell_types[i][j] == 3
                    Vbound.push_back(RowVector3d(min_x + j * cell_step_x + cell_step_x / 2,
                        min_y + i * cell_step_y + cell_step_y / 2, 0));
                }
                if (cell_types[i][j] == 2) { // for outside ones
                    Vout.push_back(RowVector3d(min_x + j * cell_step_x + cell_step_x / 2,
                        min_y + i * cell_step_y + cell_step_y / 2, 0));
                }
            }
        }
        MatrixXd VinsM = MatrixXd::Zero(Vins.size(), 3);
        MatrixXd VboundM = MatrixXd::Zero(Vbound.size(), 3);
        MatrixXd VoutM = MatrixXd::Zero(Vout.size(), 3);

        for (int i = 0; i < Vins.size();i++) {
            VinsM.row(i) = Vins[i];
        }
        for (int i = 0; i < Vout.size();i++) {
            VoutM.row(i) = Vout[i];
        }
        for (int i = 0; i < Vbound.size();i++) {
            VboundM.row(i) = Vbound[i];
        }
        viewer.data().add_points(VboundM, Eigen::RowVector3d(1, 0.5, 0));
        viewer.data().add_points(VinsM, Eigen::RowVector3d(0, 0.5, 1));
        //viewer.data().add_points(VoutM, Eigen::RowVector3d(0, 1, 0));
    }
    void draw_gridIntensity(igl::opengl::glfw::Viewer& viewer, int pickedV) {
        std::vector<RowVector3d> Vins, Vout, Vbound, C1, C2;
        for (int i = 0; i < side_length; i++) {
            for (int j = 0; j < side_length; j++) {

                if (cell_types[i][j] == 0) { // for inside ones
                    Vins.push_back(RowVector3d(min_x + j * cell_step_x,
                        min_y + i * cell_step_y, 0));
                    C1.push_back(RowVector3d(0, 0, cell_vectors[i][j][pickedV]));//cell_vectors[i][j][pickedV]));

                }
                if (cell_types[i][j] == 1) { // for boundary ones HERE
                    Vbound.push_back(RowVector3d(min_x + j * cell_step_x,
                        min_y + i * cell_step_y, 0));
                    C2.push_back(RowVector3d(0, cell_vectors[i][j][pickedV], 0));//cell_vectors[i][j][pickedV]));

                }
            }
        }
        MatrixXd VinsM = MatrixXd::Zero(Vins.size(), 3);
        MatrixXd VboundM = MatrixXd::Zero(Vbound.size(), 3);
        MatrixXd C1M = MatrixXd::Zero(Vins.size(), 3);
        MatrixXd C2M = MatrixXd::Zero(Vbound.size(), 3);


        for (int i = 0; i < Vins.size();i++) {
            VinsM.row(i) = Vins[i];
            C1M.row(i) = C1[i];
        }
        for (int i = 0; i < Vbound.size();i++) {
            VboundM.row(i) = Vbound[i];
            C2M.row(i) = C2[i];
        }
        viewer.data().add_points(VboundM, C2M);
        viewer.data().add_points(VinsM, C1M);
    }
};

void draw_curve(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V) {
    for (unsigned i = 0; i < V.rows() - 1; ++i) {
        if (V.row(i + 1) == RowVector3d(0, 0, 0)) {
            viewer.data().add_edges(V.row(i), V.row(0), Eigen::RowVector3d(0, 0, 1));
            break;
        }
        viewer.data().add_edges(V.row(i), V.row(i + 1), Eigen::RowVector3d(0, 0, 1));
    }
}
void draw_curveC(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V) {
    for (unsigned i = 0; i < V.rows(); ++i) {
        if (i == V.rows() - 1) {
            viewer.data().add_edges(V.row(i), V.row(0), Eigen::RowVector3d(0, 0, 1));
            break;
        }
        viewer.data().add_edges(V.row(i), V.row(i + 1), Eigen::RowVector3d(0, 0, 1));
    }
}

void draw_points(igl::opengl::glfw::Viewer& viewer, const MatrixXd& V, int i, int j) {
    viewer.append_mesh();
    if (i == 1) {
        viewer.data().add_points(V, Eigen::RowVector3d(1, 0, 0));
    }
    else {
        viewer.data().add_points(V, Eigen::RowVector3d(0, 1, 0));
    }
}
void draw_points2(igl::opengl::glfw::Viewer& viewer, const MatrixXd& Vm, const MatrixXd& Vc) {
    viewer.append_mesh();
    viewer.data().add_points(Vm, Eigen::RowVector3d(1, 0, 0));

    viewer.data().add_points(Vc, Eigen::RowVector3d(0, 1, 0));
}


int main(int argc, char* argv[]) {
    std::cout << "To choose vertices of mesh press 1 " << std::endl;
    std::cout << "To choose vertices of cage press 2 " << std::endl;
    std::cout << "To move vertices of cage press 3 " << std::endl;
    std::cout << "To show bound lines of mesh and cage press 4 " << std::endl;
    std::cout << "When you are finished with choosing vertices, press 5 to compute harmonic coordinates" << std::endl;
    std::cout << "To show interior and boundary cells press 6" << std::endl;

    int width = 1280;
    int height = 800;

    igl::opengl::glfw::Viewer viewer;
    Grid grid(width, height);


    viewer.core().is_animating = true;
    using RT = igl::opengl::ViewerCore::RotationType;
    viewer.core().rotation_type = RT::ROTATION_TYPE_NO_ROTATION;
    viewer.core().lighting_factor = 0;
    viewer.core().camera_zoom = 3.3699;

    MatrixXi F1;
    MatrixXd VforAligning = MatrixXd::Zero(4, 3);

    VforAligning.row(0) = RowVector3d(0, 0, 0);
    VforAligning.row(1) = RowVector3d(width, 0, 0);
    VforAligning.row(2) = RowVector3d(width, height, 0);
    VforAligning.row(3) = RowVector3d(0, height, 0);

    viewer.core().align_camera_center(VforAligning, F1);

    F1 = MatrixXi::Zero(2, 3);
    F1(0, 0) = 0; F1(0, 1) = 1; F1(0, 2) = 2;
    F1(1, 0) = 2; F1(1, 1) = 3; F1(1, 2) = 0;
    int n = 100;
    MatrixXd Vm = MatrixXd::Zero(n, 3);
    MatrixXd Vc = MatrixXd::Zero(0, 3);

    int pickedV = -1;
    bool pickingVerticesCage = false;
    bool movingVerticesCage = false;
    bool pickingVerticesMesh = false;
    bool drawCage = false;
    bool drawGrid = false;
    bool drawGridIntensity = false;
    viewer.callback_key_down =
        [&Vm, &Vc, &n, &drawGrid, &pickedV, &drawGridIntensity, &pickingVerticesCage, &pickingVerticesMesh, &movingVerticesCage, &drawCage, &VforAligning, &F1, &grid](igl::opengl::glfw::Viewer& viewer, unsigned char key, int)->bool {
        std::cout << "Key pressed: " << key << std::endl;
        if (key == '1') {
            pickingVerticesCage = false;
            movingVerticesCage = false;
            pickingVerticesMesh = true;
        }
        if (key == '2') {
            movingVerticesCage = false;
            pickingVerticesMesh = false;
            pickingVerticesCage = true;
        }
        if (key == '3') {
            pickingVerticesMesh = false;
            pickingVerticesCage = false;
            movingVerticesCage = true;
        }
        if (key == '4') {
            pickingVerticesMesh = false;
            pickingVerticesCage = false;
            movingVerticesCage = false;
            drawCage = !drawCage;
        }
        if (key == '5') {
            grid.updateTypes(Vc);
            grid.updateCoordinates();
            std::cout << "Finished computing harmonic coordinates" << std::endl;
            drawGridIntensity = !drawGridIntensity;

        }
        if (key == '6') {
            drawGrid = !drawGrid;
        }
        return false;
    };
    viewer.callback_mouse_down =
        [&Vm, &Vc, &n, &pickedV, &drawGrid, &pickingVerticesCage, &pickingVerticesMesh, &movingVerticesCage, &VforAligning, &F1](igl::opengl::glfw::Viewer& viewer, int, int)->bool {
        double x = viewer.current_mouse_x;
        double y = (viewer.core().viewport(3) - viewer.current_mouse_y);

        int fid;
        Eigen::Vector3f bc;
        if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
            viewer.core().proj, viewer.core().viewport, VforAligning, F1, fid, bc)) {
            for (int i = 0; i < n; i++) {
                if (movingVerticesCage) {
                    if (i < Vc.rows() && Vc(i, 0) >(x - 10) && Vc(i, 0) < (x + 10)
                        && Vc(i, 1) > (y - 10) && Vc(i, 1) < (y + 10)) {
                        pickedV = i;
                        break;
                    }
                    if (i == n - 1)
                        pickedV = -1;
                }
                if (!movingVerticesCage && pickingVerticesMesh) {
                    if (Vm.row(i) == RowVector3d(0, 0, 0)) {
                        Vm.row(i) = RowVector3d(x, y, 0);
                        break;
                    }
                    if (i == n - 1) {
                        std::cout << "No more points needed! " << std::endl;
                    }
                }
                if (!movingVerticesCage && pickingVerticesCage) {

                    if (i == Vc.rows()) {
                        Vc.conservativeResize(Vc.rows() + 1, Vc.cols());
                        Vc.row(i) = RowVector3d(x, y, 0);
                        break;
                    }
                    if (i == n - 1) {
                        std::cout << "No more points needed! " << std::endl;
                    }

                }
            }
            return true;
        }

        return false;
    };

    viewer.callback_mouse_up =
        [&grid, &Vm, &Vc, &n, &pickedV, &movingVerticesCage](igl::opengl::glfw::Viewer& viewer, int, int)->bool {
        double x = viewer.current_mouse_x;
        double y = (viewer.core().viewport(3) - viewer.current_mouse_y);
        if (movingVerticesCage) {
            if (pickedV != -1) {
                Vc.row(pickedV) = RowVector3d(x, y, 0);
                grid.updateMesh(Vm, Vc);

                grid.updateTypes(Vc);
                grid.updateCoordinates();

            }
        }

        return false;
    };
    viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer&)->bool {
        viewer.data().clear();
        draw_points2(viewer, Vm, Vc);
        if (drawCage) {
            draw_curveC(viewer, Vc);
            draw_curve(viewer, Vm);
        }
        if (drawGrid) {
            grid.draw_grid(viewer);
        }
        if (drawGridIntensity && pickedV != -1) {
            grid.draw_gridIntensity(viewer, pickedV);
        }
        return false;
    };
    viewer.launch(true, false, "Harmonic coordinates", width, height);

    return 0;
}