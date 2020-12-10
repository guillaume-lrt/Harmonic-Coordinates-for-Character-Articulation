#include <igl/opengl/glfw/Viewer.h>
#include <math.h>
#include <iostream>
#include <ostream>

using namespace Eigen;
using namespace std;

#define PI 3.14159265


class MeanValueCoordinates {
public:
	MeanValueCoordinates(MatrixXd& Vertices) {
		V = Vertices;
		n = Vertices.rows();
	}
	double angle(RowVector2d& P1, RowVector2d& P2, RowVector2d& P3, bool degree = false) {
		// return signed angle between P1, P2 and P3, in radian

		RowVector2d v1 = P1 - P2;
		RowVector2d v2 = P3 - P2;
		//double res = atan2(v1[0] * v2[1] - v1[1] * v2[0], v1[0] * v2[0] - v1[1] * v2[1]);
		double res = atan2(v2[1], v2[0]) - atan2(v1[1], v1[0]);
		if (res > PI) { res -= 2 * PI; }
		else if (res <= -PI) { res += 2 * PI; }
		return !degree ? res : res * 180. / PI;
	}

	double ith_weight(int i, RowVector2d &v) {
		RowVector2d v_im1 = i > 0 ? V.row(i - 1) : V.row(n - 1);
		RowVector2d v_i = V.row(i);
		//cout << V.row(i) << endl;
		RowVector2d v_ip1 = i < n - 1 ? V.row(i + 1) : V.row(0);
		double a_im1 = angle(v_im1, v, v_i);
		double a_i = angle(v_i, v, v_ip1);
		//cout << v_im1 << ", " << v << ", " << v_i << endl;
		//cout << a_im1 * 180. / PI << endl;
		double w_i = (tan(a_im1 / 2.) + tan(a_i / 2.)) / (v_i - v).norm();
		//double w_i = 1 / (v_i - v).norm();
		return w_i;
	}

	void compute_weight(RowVector2d& v) {
		// output: list of weight w.r.t v
		// w_i(v) = (tan(a_{i-1}/2) + tan(a_i/2))/norm(v_i-v) for all i
		// a_i is the angle between v_i, v and v_i+1

		weight = {};
		weight_sum = 0;

		for (int i = 0; i < n; i++) {
			double w_i = ith_weight(i,v);
			weight.push_back(w_i);
			weight_sum += w_i;
		}
	}

	void compute_lambda(RowVector2d& v) {
		// lambda_i(v) = w_i(v) / sum(w_i(v)) 
		lambda = {};
		compute_weight(v);
		for (auto& i : weight) {
			lambda.push_back(i / weight_sum);
		}
	}

	vector<double> get_lambda() {
		return lambda;
	}

private:
	MatrixXd V;
	double weight_sum;
	vector<double> weight;
	vector<double> lambda;
	int n;
};
