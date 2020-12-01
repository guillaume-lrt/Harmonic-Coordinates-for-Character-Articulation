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

enum cell_Type{UNTYPED, EXTERIOR, BOUNDARY, INTERIOR};




string to_enum(int label) {
    if(label == UNTYPED) {
        return "UNTYPED";
    }
    if(label == EXTERIOR) {
        return "EXTERIOR";
    }
    if(label == BOUNDARY) {
        return "BOUNDARY";
    }
    if(label == INTERIOR) {
        return "INTERIOR";
    }
}
struct Cell{
    vector<double> harmonicCoordinates ;
    int label;
    Cell() {
        label = UNTYPED;
    }
    void initialize(int n) {
        for (int i = 0; i<n ; i++ ) {
            harmonicCoordinates.push_back(0); //initial values are 0
        }
        label = UNTYPED; //initial label is UNTYPED
    }

    void to_string() {
        cout<<"\nCell Type: "<<to_enum(label)<<endl;
        cout<<"Harmonic Coordinates: [";

        for(auto it = harmonicCoordinates.begin(); it<harmonicCoordinates.end(); it++) {
            cout<< *it;
            if(it != harmonicCoordinates.end()-1) {
                cout<<", ";
            }
        }
        cout<<"]"<<endl;
    }
};


typedef vector<Cell> v2d;
typedef vector<v2d> Grid;



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
		-7.8, 3,
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
		7.8, 3,
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
	Model = MatrixXd(32, 2);//32

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


void createEdges(const MatrixXd& V, MatrixXd &W) {
	// shift V by 1 to draw edges
	int n = V.rows();
	for (int i = 0; i < n-1; i++) {
		W.row(i) = V.row(i + 1);
	}
	W.row(n - 1) = V.row(0);
}


void createGridCell(MatrixXd& grid, int s) {
    int squareSideCells = (int)(sqrt(pow(2,s)));

    int verticesPerSide = squareSideCells+1;
    int nbOfVertices = verticesPerSide*verticesPerSide;
    cout<< "\nsquare side nbVertices: " << squareSideCells +1;
    cout<< "\nnbVertices: " << nbOfVertices;
    int offsetX = -8;
    int offsetY = 8.5;

    int counterX = 0;
    int counterY = 0;
    int h = 2.85; //jumps between vertices
    grid = MatrixXd(nbOfVertices, 2);

    int i =0;
    while(i<nbOfVertices) {

        for(int j = 0; j<verticesPerSide; j++) {
            grid.row(i) = RowVector2d(counterX + offsetX,counterY + offsetY);
            counterX += h;
            //cout<<"\n"<<i;
            i++;

        }
        counterX = 0;
        counterY -= h;

    }

}



// ------------ main program ----------------
int main(int argc, char *argv[])
{

    igl::opengl::glfw::Viewer viewer; // create the 3d viewer


	createCage(V,E);

	MatrixXd C = MatrixXd::Zero(1, 3);			// for the color (black)
	MatrixXd W = MatrixXd::Zero(V.rows(), 2);		// V shifted by 1 to draw edges

	MatrixXd Model;
	createModel(Model);

	MatrixXd Wm = MatrixXd::Zero(Model.rows(), 2);		// edges for the model

	createEdges(V, W); //CAGE
	viewer.data().add_points(V, RowVector3d(255, 0, 0));
	viewer.data().add_edges(V, W, RowVector3d(0, 255, 0));



	createEdges(Model, Wm); // MODEL
	viewer.data().add_points(Model, RowVector3d(255, 255, 255));
	viewer.data().add_edges(Model, Wm, RowVector3d(255, 0, 255));

    MatrixXd GridVertices; //build grid
    int s = 6; //based on paper, s controls size of grid 2^s
    createGridCell(GridVertices, s);
    viewer.data().add_points(GridVertices, RowVector3d(50, 50, 0));


    int numberOrRowCells = (int)(sqrt(pow(2,s))) + 1;
    int numberOrColCells = numberOrRowCells;  //square grid we are working on
    int cageVerticesCount = V.rows();

    // Grid is a vector or vector <Cell> where Cell contains a vector of harmonic coordinates and a TYPE
    Grid grid(numberOrRowCells, v2d(numberOrColCells));

    //initialize grid cells;

    for(int i = 0; i<numberOrRowCells; i++ ){
        for(int j = 0; j<numberOrColCells; j++ ){
            grid[i][j].initialize(cageVerticesCount);
            grid[i][j].to_string(); //to print the cell type and harmonic coord
        }
    }


    viewer.data().point_size = 13; //SIZE or vertices in the viewer (circle size)
	//viewer.core(0).align_camera_center(V, F);
	viewer.launch(); // run the viewer
}
