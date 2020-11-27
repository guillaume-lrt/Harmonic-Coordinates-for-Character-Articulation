#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>
#include <math.h>

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

#define PI 3.14159265

MatrixXd V;
MatrixXi E;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	return true;
}

/**
 * Create a triangle mesh corresponding to an octagon inscribed in the unit circle
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

double angle(RowVector2d& P1, RowVector2d& P2, RowVector2d& P3, bool degree = false) {
	// return signed angle between P1, P2 and P3, in radian

	RowVector2d v1 = P1 - P2;
	RowVector2d v2 = P3 - P2;
	double res = atan2(v1[0] * v2[1] - v1[1] * v2[0], v1[0] * v2[0] - v1[1] * v2[1]);
	return !degree ? res : res * 180. / PI;
}

vector<double> compute_weight(RowVector2d &v, MatrixXd& Vertices, double &weight_sum) {
	// output: list of weight w.r.t v
	// w_i(v) = (tan(a_{i-1}/2) + tan(a_i/2))/norm(v_i-v) for all i
	// a_i is the angle between v_i, v and v_i+1

	int n = Vertices.rows();
	vector<double> weights;
	weight_sum = 0;

	for (int i = 0; i < n; i++) {
		RowVector2d v_im1 = i > 0 ? Vertices.row(i - 1) : Vertices.row(n - 1);
		RowVector2d v_i = Vertices.row(i);
		RowVector2d v_ip1 = i < n - 1 ? Vertices.row(i + 1) : Vertices.row(0);
		double a_im1 = angle(v_im1, v, v_i);
		double a_i = angle(v_i, v, v_ip1);
		double w_i = (tan(a_im1 / 2.) + tan(a_i / 2.)) / (v_i - v).norm();
		weights.push_back(w_i);
		weight_sum += w_i;
	}	
	return weights;
}

vector<double> compute_lambda(vector<double> &weight, double weight_sum) {
	// lambda_i(v) = w_i(v) / sum(w_i(v)) 
	vector<double> lambdas;
	for (auto& i : weight) {
		lambdas.push_back(i / weight_sum);
	}
	return lambdas;
}

// ------------ main program ----------------
int main(int argc, char *argv[])
{

	createCage(V,E);
	cout << V << endl;
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer

	MatrixXd C = MatrixXd::Zero(1, 3);
	MatrixXd W = MatrixXd(19, 2);		// V shifted by 1 to draw edges

	W << -2.2, 5,
		-8, 3,
		-7.2, 1,
		-2.5, 2.3,
		- 3, -4,
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
		2, 7,
		-2, 7;

	RowVector2d P1(1., 0.);
	RowVector2d P2(0., 0.);
	RowVector2d P3(-0.1, -1.);

	double a = angle(P1, P2, P3,true);

	cout << a << endl;

	viewer.data().add_points(V,C);
	viewer.data().add_edges(V, W, C);

	//viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer
}
