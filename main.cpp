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

void createModel(MatrixXd& Model) {
	Model = MatrixXd(32, 2);

	Model << 0, 6.5,
		-1, 6.3,
		-1.4, 6,
		-1.3, 5,
		-1, 4.5,
		-6.5, 2.7,
		-7, 2,
		-6, 1.9,
		-2.5, 3,
		-1.8, 2.65,
		-1.9, -1,
		-2.8, -6,
		-2.6, -6.5,
		-1.9, -6.5,
		-1.5, -6,
		-1, -3.3,
		0, -1.5,
		1, -3.3,
		1.5, -6,
		1.9, -6.5,
		2.6, -6.5,
		2.8, -6,
		1.9, -1,
		1.8, 2.65,
		2.5, 3,
		6, 1.9,
		7, 2,
		6.5, 2.7,
		1, 4.5,
		1.3, 5,
		1.4, 6,
		1, 6.3;
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

	MatrixXd C = MatrixXd::Zero(1, 3);			// for the color (black)
	MatrixXd W = MatrixXd::Zero(V.rows(), 2);		// V shifted by 1 to draw edges
	MatrixXd Vp = V;
	Vp.row(6) = RowVector2d(-4, -8);
	Vp.row(7) = RowVector2d(-1.5, -8);

	MatrixXd Wp = MatrixXd::Zero(Vp.rows(), 2);

	MatrixXd Model;
	createModel(Model);
	MatrixXd Wm = MatrixXd::Zero(Model.rows(), 2);		// edges for the model

	//RowVector2d P1(1., 0.);
	//RowVector2d P2(0., 0.);
	//RowVector2d P3(-0.1, -1.);

	MeanValueCoordinates mvc = MeanValueCoordinates(V);

	double e = 0.5;

	for (int x = -10; x <= 10; x++) {
		for (int y = -10; y <= 10; y++) {

			//int i = 7;
			double xx = x+e; double yy = y + e;
			RowVector2d P(xx, yy);
			mvc.compute_lambda(P);
			auto l = mvc.get_lambda();
			RowVector2d temp(0., 0.);
			for (int i = 0; i < l.size(); i++) {
				temp += l[i] * Vp.row(i);
			}

			//double w = mvc.ith_weight(i, P);			// compute the weight from i-th vertex of the cage to point P
			//RowVector3d T(xx, yy, w);
			viewer.data().add_edges(P, temp, C);
			viewer.data().add_points(temp,C);

		}
	}

	//cout << W << endl;
	create_edges(V, W);
	cout << W << endl;
	viewer.data().add_points(V, RowVector3d(265, 0, 0));
	viewer.data().add_edges(V, W, C);
	viewer.data().point_size = 10;

	create_edges(Vp, Wp);
	viewer.data().add_points(Vp, RowVector3d(0, 265, 265));
	viewer.data().add_edges(Vp, Wp, C);

	create_edges(Model, Wm);
	viewer.data().add_points(Model, RowVector3d(0, 0, 265));
	viewer.data().add_edges(Model, Wm, RowVector3d(0, 0, 265));

	//viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer
}
