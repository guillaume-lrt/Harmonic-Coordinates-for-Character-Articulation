#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <igl/readOFF.h>
#include <igl/writeOFF.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/jet.h>

using namespace Eigen; // to use the classes provided by Eigen library
using namespace std;

MatrixXd V;
//MatrixXi F;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
	return true;
}

/**
 * Create a triangle mesh corresponding to an octagon inscribed in the unit circle
 */
void createCage(MatrixXd &Vertices)
{
	Vertices = MatrixXd(19, 3);
	//Faces = MatrixXi(8, 3);

	Vertices << -2, 7, 0,
		-2.2, 5, 0,
		-8, 3, 0,
		-7.2, 1, 0,
		-2.5, 2.3, 0,
		-3, -4, 0,
		-3.3, -7, 0,
		-1, -7, 0,
		-0.3, -4, 0,
		0, -2.5, 0,
		0.3, -4, 0,
		1, -7, 0,
		3.3, -7, 0,
		3, -4, 0,
		2.5, 2.3, 0,
		7.2, 1, 0,
		8, 3, 0,
		2.2, 5, 0,
		2, 7, 0;
}

// ------------ main program ----------------
int main(int argc, char *argv[])
{

	createCage(V);
	cout << V << endl;
	igl::opengl::glfw::Viewer viewer; // create the 3d viewer

	MatrixXd C = MatrixXd::Zero(1, 3);
	MatrixXd W = MatrixXd(19, 3);		// V shifted by 1 to draw edges

	W << -2.2, 5, 0,
		-8, 3, 0,
		-7.2, 1, 0,
		-2.5, 2.3, 0,
		-3, -4, 0,
		-3.3, -7, 0,
		-1, -7, 0,
		-0.3, -4, 0,
		0, -2.5, 0,
		0.3, -4, 0,
		1, -7, 0,
		3.3, -7, 0,
		3, -4, 0,
		2.5, 2.3, 0,
		7.2, 1, 0,
		8, 3, 0,
		2.2, 5, 0,
		2, 7, 0,
		-2, 7, 0;

	//viewer.data().set_vertices(V);
	viewer.data().add_points(V,C);
	viewer.data().add_edges(V, W, C);
	//viewer.data().set_mesh(V, F);

	//viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer
}
