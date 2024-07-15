#include "read_initial.h"
#include "static_analysis.h"

void exportMatrixToFile(const Eigen::MatrixXd& matrix, const std::string& filename);
void result(Calculate& cal, int& pid);

int main() {

	std::string filename{ "..//..//..//fem initial set.txt" };
	std::map<std::string, std::vector<double>> variables;
	if (!readTxtFile(filename, variables)) {
		return 1;
	}

	// beam length and height and thickness(width)
	double L{};
	double h{};
	double w{};

	// mesh number
	size_t Nx{};
	size_t Ny{};

	// set material
	double rho{};
	double E{};
	double v{};

	// set external force with[apply time] and [node - ID] and [vector]
	int pid{};
	double P[2]{};

	// Boundary ID
	std::vector<size_t>	bid_list;

	// apply data
	apply_parameter(variables, L, h, w, Nx, Ny, rho, E, v, pid, P, bid_list);

	//// input material
	Material elastic{ rho, E, v };
	elastic.cross_section(h, w);

	// creat Mesh
	Mesh msh{ L, h, Nx, Ny };

	// start static analysis
	Calculate cal(msh);
	cal.apply_force_info(pid, P);
	cal.single_elem_matrix(elastic);
	cal.apply_node_f_ext(0.0);
	cal.set_M_K_F();

	// static analysis solver
	Static_Solver sol{ cal, bid_list };
	result(cal, pid);

	return 0;
}

void exportMatrixToFile(const Eigen::MatrixXd& matrix, const std::string& filename) {

	std::ofstream file(filename);
	if (file.is_open()) {
		for (int i = 0; i < matrix.rows(); ++i) {
			for (int j = 0; j < matrix.cols(); ++j) {
				file << matrix(i, j);

				if (j < matrix.cols() - 1) {
					file << "\t";
				}
			}

			file << std::endl;
		}

		file.close();
		std::cout << "Matrix has been exported to " << filename << std::endl;
	}

	else {
		std::cerr << "Unable to open file " << filename << " for writing!" << std::endl;
	}
}

void result(Calculate& cal, int& pid) {

	std::cout << "\nnumerical solution = " << cal.msh.nodes[pid].displacement(1) << std::endl;
	std::cin.get();

}