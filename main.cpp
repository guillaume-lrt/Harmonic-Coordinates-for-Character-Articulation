#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include <ostream>
#include <math.h>

using namespace Eigen;
using namespace std;
enum cell_Type{UNTYPED, EXTERIOR, BOUNDARY, INTERIOR};


//const double h = 2; //jumps between cells of the cage 2.85
double h = 1; //jumps between cells of the cage 2.85
const int interpolationPrecision = 50; //10 point on edges to detect cells
int s = 10; //based on paper, s controls size of grid 2^s
int vertex_count = 5;   // vertex to move with keys 1,2,3 and 4
                        // use 'p' and 'm' to change vertex


igl::opengl::glfw::Viewer viewer; // create the 3d viewer

MatrixXd Cage; //Cage vertices
MatrixXi CageEdgesIndices; //Cage Edges indices!!
MatrixXd Ec; //Cage edges

MatrixXd init_Model;
MatrixXd Model; //our ineterior little 2D human
MatrixXd Em; //Model edges

MatrixXd GridVertices; //Grid 2^s cells behind the model and cage

double offsetX = -8; //X offset of the cage -8
//const double offsetY = 8; //Y offset of the cage to enclose model 8.5
double offsetY = 9; //Y offset of the cage to enclose model 8.5



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
    return "Not Found";
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

Grid grid;

void printVector(vector<double> v) {
	cout << "[";
	for (auto& i : v) {
		cout << i << "  ";
	}
	cout << "]" << endl;
}


/**
 * Create a cage
 */
void createHumanCage(MatrixXd &Vertices, MatrixXi &Edges)
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

void createHumanModel(MatrixXd& Model) {
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

void createBiggerHumanCage(MatrixXd& Vertices, MatrixXi& Edges)
{

    Vertices = MatrixXd(19, 2);
    Edges = MatrixXi(19, 2);

    Vertices << -8., 28., -8.8, 20., -31.2, 12., -28.8, 4., -10.,
        9.2, -12., -16., -13.2, -28., -4., -28., -1.2, -16.,
        0., -10., 1.2, -16., 4., -28., 13.2, -28., 12.,
        -16., 10., 9.2, 28.8, 4., 31.2, 12., 8.8, 20.,
        8., 28.;

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

void createBiggerHumanCage_pos1(MatrixXd& Vertices)
{

    Vertices = MatrixXd(19, 2);

    Vertices << -8., 28., 
        -8.8, 20., 
        -31.2, 12., 
        -28.8, 4., 
        -10., 9.2, 
        -16., -16., 
        -19.2, -32.,
        -9, -32., 
        -3.2, -16.,
        0., -10., 
        1.2, -16., 
        4., -28., 
        13.2, -28., 
        12.,-16., 
        10., 9.2, 
        28.8, 4., 
        31.2, 12., 
        8.8, 20.,
        8., 28.;
}

void createBiggerHumanCage_pos2(MatrixXd& Vertices)
{

    Vertices = MatrixXd(19, 2);

    Vertices << -8., 28.,
        -8.8, 20.,
        -31.2, 12.,
        -28.8, 4.,
        -10., 9.2,
        -16., -16.,
        -17.2, -21.,
        -5, -21.,
        -3.2, -16.,
        0., -10.,
        1.2, -16.,
        4., -28.,
        13.2, -28.,
        12., -16.,
        10., 9.2,
        28.8, 4.,
        31.2, 12.,
        8.8, 20.,
        8., 28.;
}

void createBiggerHumanCage_head(MatrixXd& Vertices)
{

    Vertices = MatrixXd(19, 2);

    Vertices << -12., 32.,
        -8, 19.,
        -31.2, 12.,
        -28.8, 4.,
        -10., 9.2,
        -12., -16., 
        -13.2, -28., 
        -4., -28.,
        -1.2, -16.,
        0., -10.,
        1.2, -16.,
        4., -28.,
        13.2, -28.,
        12., -16.,
        10., 9.2,
        28.8, 4.,
        31.2, 12.,
        8.8, 20.,
        4., 32.;
}

void createBiggerHumanModel(MatrixXd& Model) {
    Model = MatrixXd(32, 2);//32

    Model << 0., 26., -4., 25.2, -5.6, 24., -5.2, 20., -4.,
        18., -26., 10.8, -28., 8., -24., 7.6, -10., 12.,
        -7.2, 10.6, -7.6, -4., -11.2, -24., -10.4, -26., -7.6,
        -26., -6., -24., -4., -13.2, 0., -6., 4., -13.2,
        6., -24., 7.6, -26., 10.4, -26., 11.2, -24., 7.6,
        -4., 7.2, 10.6, 10., 12., 24., 7.6, 28., 8.,
        26., 10.8, 4., 18., 5.2, 20., 5.6, 24., 4.,
        25.2;
}

void createSquareCage(MatrixXd& Vertices, MatrixXi& Edges) {
    Vertices = MatrixXd(8, 2);

    Vertices << -3.5, -3.5,
        -3.5, -0.5,
        -3.5, 2.5,
        -0.5, 2.5,
        2.5, 2.5,
        2.5, -0.5,
        2.5, -3.5,
        -0.5, -3.5;

    Edges = MatrixXi(8, 2);

    Edges << 0, 1,
        1, 2,
        2, 3,
        3, 4,
        4, 5,
        5, 6,
        6, 7,
        7, 0;
}

void createSquareModel(MatrixXd& Model) {
    Model = MatrixXd(8, 2);

    Model << -1.75, -1.75,
        -1.75, -0.25,
        -1.75, 1.25,
        -0.25, 1.25,
        1.25, 1.25,
        1.25, -0.25,
        1.25, -1.75,
        -0.25, -1.75;
}

void createSquareCageS10(MatrixXd& Vertices, MatrixXi& Edges) {
    Vertices = MatrixXd(8, 2);

    Vertices << -6.5, -17.5,
        -5.5, -4.5,
        -4.5, 8.5,
        8.5, 7.5,
        21.5, 6.5,
        20.5, -6.5,
        19.5, -19.5,
        6.5, -18.5;

    Edges = MatrixXi(8, 2);

    Edges << 0, 1,
        1, 2,
        2, 3,
        3, 4,
        4, 5,
        5, 6,
        6, 7,
        7, 0;
}

void createSquareModelS10(MatrixXd& Model) {
    Model = MatrixXd(8, 2);

    Model << 2.25, -12.75, 2.25, -8.25, 2.25, -3.75, 6.75, -3.75, 11.25, -3.75, 11.25, -8.25, 11.25, -12.75, 6.75, -12.75;
}

void createSquareCageS12(MatrixXd& Vertices, MatrixXi& Edges) {
    Vertices = MatrixXd(8, 2);

    Vertices << -6.5, -45.5,
        -5.5, -4.5,
        -4.5, 8.5,
        8.5, 7.5,
        50.5, 6.5,
        49.5, -9.5,
        48.5, -47.5,
        6.5, -46.5;

    Edges = MatrixXi(8, 2);

    Edges << 0, 1,
        1, 2,
        2, 3,
        3, 4,
        4, 5,
        5, 6,
        6, 7,
        7, 0;
}

void createSquareModelS12(MatrixXd& Model) {
    Model = MatrixXd(8, 2);

    Model << 2.25, -12.75, 2.25, -8.25, 2.25, -3.75, 6.75, -3.75, 11.25, -3.75, 11.25, -8.25, 11.25, -12.75, 6.75, -12.75;
}

void createEdges(const MatrixXd& V, MatrixXd &W) {
	// shift V by 1 to draw edges
	int n = V.rows();
	for (int i = 0; i < n-1; i++) {
		W.row(i) = V.row(i + 1);
	}
	W.row(n - 1) = V.row(0);
}


void createGridVertices(MatrixXd& grid, int s) {
    int squareSideCells = (int)(sqrt(pow(2,s)));

    int verticesPerSide = squareSideCells+1;
    int nbOfVertices = verticesPerSide*verticesPerSide;
    cout<< "\nsquare side nbVertices: " << squareSideCells +1;
    cout<< "\nnbVertices: " << nbOfVertices;

    double counterX = 0;
    double counterY = 0;

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

    int x = 0; //x in terms of cell --> remember x goes vertically down as row and y horiz as columns
    int y = 0; //y in terms of cell
    for(int i =0; i<Cage.rows(); i++) {
        mapVerticesToGridCoord( x, y, i);
        cout<<"\n x"<<x<<" y"<<y<<" i"<<i;
        //mark the x,y cell as boundary
        if(grid[x][y].label != BOUNDARY) {
            grid[x][y].label = BOUNDARY;
            grid[x][y].harmonicCoordinates[i] = 1; //harmonic coordinates of this cell is 0...1..0 where 1 is at ith index
            grid[x][y].to_string();
            //viewer.data().add_label(GridVertices.row(9*x+y),"BOUNDARY");
        }
        //viewer.data().add_label(Cage.row(i),"BOUNDARY");
    }

    //Loop over edges of cage, interpolate with certain precision and mark grid boundaries

    int v1Index = 0;
    int v2Index = 0;

    double alpha = 1.0/interpolationPrecision;
    VectorXd v;

    for(int i=0; i<CageEdgesIndices.rows(); i++) {
        v1Index = CageEdgesIndices(i,0);
        v2Index = CageEdgesIndices(i,1);
        //cout<<endl<<"Interpolation between: "<<v1Index<<" "<<v2Index<<endl;
        double interpV1 = 0;
        double interpV2 = 0;

        for(int j = 0; j<interpolationPrecision; j++) {
            // the next two interpolation values should be used as harmonic coordinates for these boundary cells
            interpV1 = j*alpha;
            interpV2 = (1-interpV1);
            v = interpV1*Cage.row(v1Index) + interpV2*Cage.row(v2Index); //linear interpolations

            mapVerticesToGridCoord(x, y, v(0), v(1));

            if(grid[x][y].label != BOUNDARY) {
                grid[x][y].label = BOUNDARY;
                grid[x][y].harmonicCoordinates[v1Index] = interpV1;
                grid[x][y].harmonicCoordinates[v2Index] = interpV2;
                grid[x][y].to_string();
                //viewer.data().add_label(GridVertices.row(9*x+y),"BOUNDARY");
            }
        }
    }
}

vector<int> grid_neighbour(Grid& grid,int i, int j) {
	// output: vector of indices of neighbour of grid[i][j]
	vector<int> res;
	int n = grid.size();

	for (int k = -1; k <= 1; k++) {
		for (int l = -1; l <= 1; l++) {
			if (k != 0 || l != 0) {
                if (k == 0 || l == 0) {
                    if (k + i <= n - 1 && k + i >= 0 && l + j <= n - 1 && l + j >= 0) {
                        res.push_back(k + i);
                        res.push_back(l + j);
                    }
                }
			}
		}
	}
	return res;
}

void exploreGrid(Grid& grid, int i, int j) {

	if (grid[i][j].label == UNTYPED) {

		grid[i][j].label = EXTERIOR;
        //viewer.data().add_label(GridVertices.row(9*i+j),"EXTERIOR");
		auto neigh = grid_neighbour(grid, i, j);
		for (int k = 0; k < neigh.size(); k = k + 2) {
            exploreGrid(grid, neigh[k], neigh[k + 1]);
		}
	}
}


void labelGrid(Grid& grid) {
    fillBoundaryCells(grid); //Fill BOUNDARY edges with interpolation and then propagate for exterior
    cout<<"\nBoundary Labeled";
    int n = grid.size();

	exploreGrid(grid, 0, 0); //TODO can be changed to something more robust
    exploreGrid(grid, 0, n - 1);
    exploreGrid(grid, n - 1, 0);
    exploreGrid(grid, n-1, n-1);

    for (auto& i : grid) { // Everything that is not labeled now can be considered interior
		for (auto& j : i) {
			if (j.label == UNTYPED) j.label = INTERIOR;
		}
	}
}

void printGrid(Grid& grid) {
    cout<<"\n\n";
    cout<<"----------------- Printing the Grid ----------------- \n\n";
    int n = grid.size();
    for(int i=0; i<n; i++) {
        for(int j=0; j<n; j++) {
            if(grid[i][j].label == UNTYPED) cout<<"U";
            else if(grid[i][j].label == EXTERIOR) cout<<"-  ";
            else if(grid[i][j].label == BOUNDARY) cout<<"|  ";
            else if(grid[i][j].label == INTERIOR) cout<<"*  ";
        }
        cout<<endl;
    }
    cout<<endl;
}

double forceZeroLaplacian(Grid &grid, int i, int j) {

    //Need to force zero laplacian following the discrete law of laplacian and return the current change from old to new

    //assume a cage is all over the model and thus we will always have surrouding 4 cells ! this can be easly changed later
    if(i<0 || j<0) return 0;

    int n = grid.size(); //assume square grid
    int harmoniceSize = grid[i][j].harmonicCoordinates.size(); //assume square grid

    double holder = 0;
    double change = 0;
    for(int v = 0; v<harmoniceSize; v++) {
        //formula is [i-1][j] + [i+1][j] + [i][j-1]  + [i][j+1] = 4  * [i][j]
        holder = grid[i][j].harmonicCoordinates[v];
        grid[i][j].harmonicCoordinates[v] = grid[i-1][j].harmonicCoordinates[v] + grid[i+1][j].harmonicCoordinates[v]
                                            + grid[i][j-1].harmonicCoordinates[v] + grid[i][j+1].harmonicCoordinates[v];

        grid[i][j].harmonicCoordinates[v] /= 4.0;

        change += abs(grid[i][j].harmonicCoordinates[v] - holder);
    }

    return change/harmoniceSize;
}


void propagateLaplacian(Grid &grid, double threshold) {

    int n = grid.size();
    double maxChange = 0; //positive always, change in harmonic values as mean
    double currentChange = 0; // also positive, local change of harmonic values as mean

    int counter = 0;
    while(true) { //till we reach the thresh we need
        maxChange = 0;
        counter++;
        for(int i=0; i<n; i++) {
            for (int j = 0; j < n; j++) {
                if(grid[i][j].label == INTERIOR) { //the propagation only happens on interior cells and in order/efficiency
                    currentChange = forceZeroLaplacian(grid, i, j);
                    if(currentChange > maxChange) {
                        maxChange = currentChange; //TODO later save the current change over each cell and present them
                    }
                }
            }
        }
        //cout<<"Max Change: " << maxChange << "    Round: " << counter << endl;
        if(maxChange < threshold) break;
    }

}


void updateModelBasedOnHCoordinates(MatrixXd &Model, MatrixXd  &Cage, Grid &grid) {
    //Loop over model and map to cage cell to update coordinates
    int x = 0;
    int y = 0;
    for(int i=0; i<Model.rows(); i++) {
        mapVerticesToGridCoord(x, y, Model(i,0), Model(i,1));

        for(int c=0; c<Cage.rows(); c++) {
            if(c==0) {
                Model.row(i) = grid[x][y].harmonicCoordinates[c]*Cage.row(c);
                continue;
            }
            Model.row(i) += grid[x][y].harmonicCoordinates[c]*Cage.row(c);
        }
    }
}

void refreshViewer() {

    viewer.data().clear();

    Model = init_Model;             // to avoid glitches, start from initial model
    
    updateModelBasedOnHCoordinates(Model, Cage, grid);
    //add Cage and Cage edges
    viewer.data().add_points(Cage, RowVector3d(255, 0, 0));
    createEdges(Cage, Ec);
    viewer.data().add_edges(Cage, Ec, RowVector3d(0, 255, 0));

    //add Model and model edges
    viewer.data().add_points(Model, RowVector3d(255, 255, 255));
    createEdges(Model, Em);
    viewer.data().add_edges(Model, Em, RowVector3d(255, 0, 255));

    //add te Grid layer
    viewer.data().add_points(GridVertices, RowVector3d(50, 50, 0));


    //viewer.data().point_size = 4; //SIZE or vertices in the viewer (circle size)
    viewer.data().point_size = 5;
    viewer.data().line_width = 2;
    viewer.data().show_custom_labels = true;
}



// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
    cout << "pressed Key: " << key << " " << (unsigned int)key << std::endl;
    if (key == 'F')
    {
        Cage(vertex_count,0) -= 1;
        //createBiggerHumanModel(Model);
        cout << Cage(vertex_count, 0) << endl;
        refreshViewer();
    }
    if (key == 'H')
    {
        Cage(vertex_count,0) += 1;
        //createBiggerHumanModel(Model);
        cout << Cage(vertex_count, 0) << endl;
        refreshViewer();
    }

    if (key == 'G')
    {
        Cage(vertex_count,1) -= 1;
        //createBiggerHumanModel(Model);
        cout << Cage(vertex_count, 1) << endl;
        refreshViewer();
    }
    if (key == 'T')
    {
        Cage(vertex_count,1) += 1;
        //createBiggerHumanModel(Model);
        cout << Cage(vertex_count, 1) << endl;
        refreshViewer();
    }

    if (key == 'P') {
        cout << vertex_count << endl;
        vertex_count++;
        cout << vertex_count << endl;
    }
    if (key == 'M' && vertex_count > 0) {
        vertex_count--;
    }

    if (key == '1') {
        // move left leg down
        createBiggerHumanCage_pos1(Cage);
        refreshViewer();
    }

    if (key == '2') {
        // move left leg up
        createBiggerHumanCage_pos2(Cage);
        refreshViewer();
    }

    if (key == '3') {
        // move head up left
        createBiggerHumanCage_head(Cage);
        refreshViewer();
    }

    if (key == 'R') {
        // reset the view
        createBiggerHumanCage(Cage, CageEdgesIndices);
        //createBiggerHumanModel(Model);
        
        refreshViewer();
    }
    return false;
}

// ------------ main program ----------------
int main(int argc, char *argv[]) {



    //FOR normal Human Cage uncomment this
    //h = 2;
    //s = 6;
    //offsetY = 8;
    //createHumanCage(Cage, CageEdgesIndices); //Cage Edges indices are needed Keep them !!
    //createHumanModel(Model);

    h = 1;
    s = 13;
    offsetY = 35;
    offsetX = -35;
    createBiggerHumanCage(Cage, CageEdgesIndices);
    createBiggerHumanModel(Model);


    //For square normal model
    /*h=2;
    s=6;
    offsetY = 8;
    createSquareCage(Cage, CageEdgesIndices);
    createSquareModel(Model);*/


    //For square model for s = 10
    /*h=1;
    s=10;
    offsetY = 9;
    createSquareCageS10(Cage, CageEdgesIndices);
    createSquareModelS10(Model);*/



    //For square model s = 12
    /*h=1;
    s=12;
    offsetY = 9;
    createSquareCageS12(Cage, CageEdgesIndices);
    createSquareModelS12(Model);*/

    init_Model = Model;


	Ec = MatrixXd::Zero(Cage.rows(), 2);		// CAGE shifted by 1 to draw edges
	createEdges(Cage, Ec); //CAGE


    Em = MatrixXd::Zero(Model.rows(), 2);		// edges for the model
    createEdges(Model, Em); // MODEL

	//build grid of cells
    createGridVertices(GridVertices, s); //Grid vertices are the vertices the appear in the background
    //viewer.data().add_points(GridVertices, RowVector3d(50, 50, 0));

    int numberOrRowCells = (int)(sqrt(pow(2,s)));
    int numberOrColCells = numberOrRowCells;  //square grid we are working on
    int cageVerticesCount = Cage.rows();

    // Grid is a vector or vector <Cell> where Cell contains a vector of harmonic coordinates and a TYPE
    Grid grid_(numberOrRowCells, v2d(numberOrColCells));
    grid = grid_;
    //initialize grid cells;
    for(int i = 0; i<numberOrRowCells; i++ ){
        for(int j = 0; j<numberOrColCells; j++ ){
            grid[i][j].initialize(cageVerticesCount);
            //grid[i][j].to_string(); //to print the cell type and harmonic coord
        }
    }


    //label the cage cells as INTERIOR, BOUNDARY, EXTERIOR
    labelGrid(grid);
    cout<<"\nGRID Labeled";
    //double threshold = 0.00001; // 10^-5 as paper used
    double threshold = 0.0000000001; // 10^-10
    //Propagate laplacian
    propagateLaplacian(grid, threshold);
    cout<<"\nLaplacina Propagated";
    //Update Current positions based on HC
    //updateModelBasedOnHCoordinates(Model, Cage, grid);

    //To see the grid in a nice way
    printGrid(grid);
    int interiorCellCount = 0;
    for(int i = 0; i<numberOrRowCells; i++ ){
        for(int j = 0; j<numberOrColCells; j++ ){
            // grid[i][j].to_string(); //to print the cell type and harmonic coord
            if(grid[i][j].label == INTERIOR) interiorCellCount++;
            for(auto e: grid[i][j].harmonicCoordinates) {
                if(e<0) {
                    std::cout<<"ERRORR NEGATIVE VALUE";
                    break;
                }
            }
        }
    }

    cout<<"\nINTERIOR Cell Count: "<<interiorCellCount<<"\n";
    //Refresh the viewer to show results
    refreshViewer(); //Takes care of rebuilding the edges also

    viewer.callback_key_down = &key_down;
    viewer.core().is_animating = true; // animation is active
    viewer.launch(); // run the viewer
}
