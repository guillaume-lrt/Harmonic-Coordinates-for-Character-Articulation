//
// Created by ghosto on 12/1/20.
//

#include <igl/opengl/glfw/Viewer.h>
#include <math.h>

using namespace Eigen;
using namespace std;

typedef vector<double> v1d;
typedef vector<v1d> v2d;
typedef vector<v2d> v3d;

class HarmonicCoordinates {
public:
    HarmonicCoordinates(MatrixXd& Vertices) {
        V = Vertices;
        n = Vertices.rows();
    }


private:
    MatrixXd V;
    int n;
};