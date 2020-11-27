#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <math.h>
#include "MeanValueCoordinates.cpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace Eigen;
using namespace std;

#define PI 3.14159265

MatrixXd V;
MatrixXi E;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	return true;
}

void print(vector<double> v) {
	cout << "[";
	for (auto& i : v) {
		cout << i << "  ";
	}
	cout << "]" << endl;
}

/**
 * Create a cage
 */
void createCage(MatrixXd &Vertices, MatrixXi &Edges)
{
	Vertices = MatrixXd(19, 2);
	Edges = MatrixXi(19, 2);

	Vertices << -2, 7, 
		-2.2, 5, 
		-8, 3, 
		-7.2, 1, 
		-2.5, 2.3, 
		-3, -4, 
		-3.3, -7,
		-1, -7,
		-0.3, -4,
		0, -2.5, 
		0.3, -4, 
		1, -7, 
		3.3, -7,
		3, -4,
		2.5, 2.3,
		7.2, 1, 
		8, 3,
		2.2, 5,
		2, 7;

	Edges << 0, 1,
		1, 2,
		2, 3,
		3, 4,
		4, 5,
		5, 6,
		6, 7,
		7, 8,
		8, 9,
		9, 10,
		10, 11,
		11, 12,
		12, 13,
		13, 14,
		14, 15,
		15, 16,
		16, 17,
		17, 18,
		18, 0;
}

void create_edges(const MatrixXd& V, MatrixXd &W) {
	// shift V by 1 to draw edges
	int n = V.rows();
	for (int i = 0; i < n-1; i++) {
		W.row(i) = V.row(i + 1);
	}
	W.row(n - 1) = V.row(0);
}

// ------------ main program ----------------
int main(int argc, char *argv[])
{

	createCage(V,E);
	cout << V << endl;
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer

	MatrixXd C = MatrixXd::Zero(1, 3);			// color 
	MatrixXd W = MatrixXd::Zero(V.rows(), 2);		// V shifted by 1 to draw edges

	RowVector2d P1(1., 0.);
	RowVector2d P2(0., 0.);
	RowVector2d P3(-0.1, -1.);

	MeanValueCoordinates mvc = MeanValueCoordinates(V);
	//mvc.compute_lambda(P2);
	//auto l = mvc.get_lambda();
	//print(l);

	for (int x = -8; x <= 8; x++) {
		for (int y = -7; y <= 7; y++) {
			int i = 0;
			RowVector2d P(x, y);
			double w = mvc.ith_weight(i, P);			// compute the weight from i-th vertex of the cage to point P
			RowVector3d T(x, y, w);
			viewer.data().add_points(T,C);
		}
	}

	//cout << W << endl;
	create_edges(V, W);
	cout << W << endl;
	viewer.data().add_points(V,C);
	viewer.data().add_edges(V, W, C);
	viewer.data().point_size = 10;

	//viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer
}
