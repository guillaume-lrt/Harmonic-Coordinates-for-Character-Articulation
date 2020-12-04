#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <math.h>


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

igl::opengl::glfw::Viewer viewer; // create the 3d viewer

MatrixXd Cage; //Cage vertices
MatrixXi CageEdges; //Cage Edges

MatrixXd Model; //our ineterior little 2D human
MatrixXd Em; //Model edges

MatrixXd GridVertices; //Grid 2^s cells behind the model and cage

const int offsetX = -8; //X offset of the cage
const int offsetY = 8.5; //Y offset of the cage to enclose model
const int h = 2.85; //jumps between cells of the cage

const int interpolationPrecision = 20; //10 point on edges to detect cells

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

    int counterX = 0;
    int counterY = 0;

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

vector<int> grid_neighbour(Grid& grid,int i, int j) {
	// output: vector of indices of neighbour of grid[i][j]
	vector<int> res;
	int n = grid.size();
	//cout << "n: " << n << ", " << grid[0].size() << endl;
	for (int k = -1; k <= 1; k++) {
		for (int l = -1; l <= 1; l++) {
			if (k != 0 || l != 0) {
				if (k + i <= n-1 && k + i >= 0 && l + j <= n-1 && l + j >= 0) {
					res.push_back(k + i);
					res.push_back(l + j);
				};
			}
		}
	}
	return res;
}

void explore_grid(Grid& grid, int i, int j) {
	//cout << i << ", " << j << endl;
	if (grid[i][j].label == UNTYPED) {
		//cout << "test\n";
		grid[i][j].label = EXTERIOR;
		auto neigh = grid_neighbour(grid, i, j);
		for (int k = 0; k < neigh.size(); k = k + 2) {
			explore_grid(grid, neigh[k], neigh[k + 1]);
		}
	}
}


void label_cage(Grid& grid) {
	int n = grid.size();
	explore_grid(grid, 0, 0);
	explore_grid(grid, 0, n - 1);
	explore_grid(grid, n - 1, 0);
	explore_grid(grid, n-1, n-1);
	for (auto& i : grid) {
		for (auto& j : i) {
			if (j.label == UNTYPED) j.label = INTERIOR;
		}
	}
}

void mapVerticesToGridCoord(int &x, int &y, int vCageindex) {
    y = (int)floor((Cage(vCageindex,0) - offsetX)/h); //FLIP y x for grid
    x = (int)floor((offsetY - Cage(vCageindex,1))/h); // y is opposite sign in scree cooridnates
}

void mapVerticesToGridCoord(int &xCell, int &yCell, double xVertix, double yVertix) {
    yCell = (int)floor((xVertix - offsetX)/h); //FLIP y x for grid
    xCell = (int)floor((offsetY - yVertix)/h); // y is opposite sign in scree cooridnates
}

void fillBoundaryCells(Grid &grid) {
    //Loop over cage vertices and map them to cells
    //interpolate and fill other boundary cells
    // to do the mapping to the cage graph we need to do a translation of
    // offsetX and offsetY for each point we have before filling the grid
    cout<<"Entered Boundary cells "<<grid.size()<< " "<<Cage.rows();

    int x = 0; //x in terms of cell --> remember x goes vertically down as row and y horiz as columns
    int y = 0; //y in terms of cell
    for(int i =0; i<Cage.rows(); i++) {
        mapVerticesToGridCoord( x, y, i);
        cout<<"\n"<<x<<" "<<y;
        //mark the x,y cell as boundary
        if(grid[x][y].label != BOUNDARY) {
            grid[x][y].label = BOUNDARY;
            viewer.data().add_label(GridVertices.row(9*x+y),"BOUNDARY");
        }
        //viewer.data().add_label(Cage.row(i),"BOUNDARY");
    }

    //Loop over edges of cage, interpolate with certain precision and mark grid boundaries

    int v1Index = 0;
    int v2Index = 0;

    double alpha = 1.0/interpolationPrecision;
    VectorXd v;

    for(int i=0; i<CageEdges.rows(); i++) {
        v1Index = CageEdges(i,0);
        v2Index = CageEdges(i,1);
        //cout<<endl<<"Interpolation between: "<<v1Index<<" "<<v2Index<<endl;
        for(int j = 0; j<interpolationPrecision; j++) {
            // the next two interpolation values should be used as harmonic coordinates for these boundary cells
            v = j*alpha*Cage.row(v1Index) + (1-j*alpha)*Cage.row(v2Index); //linear interpolations

            mapVerticesToGridCoord(x, y, v(0), v(1));

            if(grid[x][y].label != BOUNDARY) {
                grid[x][y].label = BOUNDARY;
                viewer.data().add_label(GridVertices.row(9*x+y),"BOUNDARY");
            }
        }
        //cout<<endl;
    }

}

// ------------ main program ----------------
int main(int argc, char *argv[])
{

	createCage(Cage, CageEdges);

	MatrixXd C = MatrixXd::Zero(1, 3);			// for the color (black)
	MatrixXd W = MatrixXd::Zero(Cage.rows(), 2);		// CAGE shifted by 1 to draw edges


	createModel(Model);

	Em = MatrixXd::Zero(Model.rows(), 2);		// edges for the model

	createEdges(Cage, W); //CAGE
	viewer.data().add_points(Cage, RowVector3d(255, 0, 0));
	viewer.data().add_edges(Cage, W, RowVector3d(0, 255, 0));



	createEdges(Model, Em); // MODEL
	viewer.data().add_points(Model, RowVector3d(255, 255, 255));
	viewer.data().add_edges(Model, Em, RowVector3d(255, 0, 255));

     //build grid
    int s = 6; //based on paper, s controls size of grid 2^s
    createGridCell(GridVertices, s);
    viewer.data().add_points(GridVertices, RowVector3d(50, 50, 0));


    int numberOrRowCells = (int)(sqrt(pow(2,s))) + 1;
    int numberOrColCells = numberOrRowCells;  //square grid we are working on
    int cageVerticesCount = Cage.rows();

    // Grid is a vector or vector <Cell> where Cell contains a vector of harmonic coordinates and a TYPE
    Grid grid(numberOrRowCells, v2d(numberOrColCells));

    //initialize grid cells;

	grid[2][0].label = BOUNDARY;
	grid[3][0].label = BOUNDARY;
	grid[4][0].label = BOUNDARY;
	grid[5][0].label = BOUNDARY;
	grid[2][1].label = BOUNDARY;
	grid[3][1].label = BOUNDARY;
	grid[4][1].label = BOUNDARY;
	grid[5][1].label = BOUNDARY;
	grid[2][2].label = BOUNDARY;
	grid[3][2].label = BOUNDARY;
	grid[4][2].label = BOUNDARY;
	grid[5][2].label = BOUNDARY;
	grid[2][3].label = BOUNDARY;

	grid[5][3].label = BOUNDARY;
	grid[2][4].label = BOUNDARY;

	grid[5][4].label = BOUNDARY;
	grid[2][5].label = BOUNDARY;

	grid[5][5].label = BOUNDARY;
	grid[2][6].label = BOUNDARY;

	grid[5][6].label = BOUNDARY;
	grid[2][7].label = BOUNDARY;
	grid[3][7].label = BOUNDARY;
	grid[4][7].label = BOUNDARY;
	grid[5][7].label = BOUNDARY;

	grid[0][4].label = BOUNDARY;
	grid[0][5].label = BOUNDARY;
	grid[0][6].label = BOUNDARY;
	grid[1][4].label = BOUNDARY;
	grid[1][5].label = BOUNDARY;
	grid[1][6].label = BOUNDARY;

	grid[7][4].label = BOUNDARY;
	grid[7][5].label = BOUNDARY;
	grid[7][6].label = BOUNDARY;
	grid[6][4].label = BOUNDARY;
	grid[6][5].label = BOUNDARY;
	grid[6][6].label = BOUNDARY;

	label_cage(grid);

    for(int i = 0; i<numberOrRowCells; i++ ){
        for(int j = 0; j<numberOrColCells; j++ ){
            //grid[i][j].initialize(cageVerticesCount);
			cout << "\nCell[" << i << "][" << j << "]";
            grid[i][j].to_string(); //to print the cell type and harmonic coord
        }
    }


    //fill the boundary edges first
    cout<<"Calling boundaries:";
    fillBoundaryCells(grid);
    cout<<"\nDONE";
    viewer.data().point_size = 13; //SIZE or vertices in the viewer (circle size)
	//viewer.core(0).align_camera_center(V, F);
    viewer.data().show_custom_labels = true;
	viewer.launch(); // run the viewer
}
