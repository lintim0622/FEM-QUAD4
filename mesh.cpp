#include "mesh.h"

// class Material
Material::Material() : _rho{ 0.0 }, E{ 0.0 }, v{ 0.0 } {

	h = 0.0;
	_w = 0.0;
	_A = 0.0;
	I = 0.0;

	E1 = 0.0;
	E2 = 0.0;
	G = 0.0;
	_D.setZero();

}

Material::Material(double rho, double E, double v) : _rho{ rho }, E{ E }, v{ v } {

	h = 0.0;
	_w = 0.0;
	_A = 0.0;
	I = 0.0;

	E1 = E / (1.0 - v * v);
	E2 = v * E / (1.0 - v * v);
	G = E / (2.0 * (1.0 + v));
	_D << E1,  E2,  0.0,
		  E2,  E1,  0.0,
		  0.0, 0.0, G;

}

Material::~Material() {}

double& Material::rho() { return _rho; }

double& Material::w() { return _w; }

double& Material::A() { return _A; }

Eigen::Matrix<double, 3, 3>& Material::D() { return _D; }

void Material::inputData() {

	// material parameter
	// rho -> density, E -> Young's moduls, v -> poisson ratio
	/*double rho = 7.964e-9; double E = 200e+9; double v = 0.31;*/
	std::cout << "\nplease material parameter:\n";
	std::cout << "density = ";
	std::cin >> _rho;
	std::cout << "Young's moduls = ";
	std::cin >> E;
	std::cout << "poisson ratio = ";
	std::cin >> v;

}

void Material::cross_section(double h, double w) {

	this->h = h;
	_w = w;
	_A = w * h;
	I = w * h * h * h / 12.0;

}

void Material::outputFile() {

	// output file
	std::ofstream outputFile("mesh and material data.txt", std::ios::out);

	if (!outputFile.is_open()) {
		std::cout << "The file open is failed\n";
	}

	outputFile << "# ---------- material ----------" << std::endl;
	std::cout.fill(' ');
	outputFile << "rho" << std::setw(15) << std::right << _rho << std::endl;
	outputFile << "E" << std::setw(17) << std::right << E << std::endl;
	outputFile << "v" << std::setw(17) << std::right << v << std::endl;
	outputFile << std::endl;

	outputFile.close();

	std::cout << "\nFinished!\n";

}


// class Node
Node::Node() {

	nid = 0;
	gid[0] = 0;
	gid[1] = 0;
	position.setZero();
	displacement.setZero();
	stress.setZero();
	strain.setZero();
	f_ext.setZero();

}

Node::~Node() {}


// class Element
Element::Element() : n1{nullptr}, n2{ nullptr }, n3{ nullptr }, n4{ nullptr } {

	eid = 0;
	me.setZero();
	ke.setZero();

}

Element::~Element() {}


// class Mesh
Mesh::Mesh(double L, double h, size_t Nx, size_t Ny) {

	size_t Nxl = Nx + 1;
	size_t Nyl = Ny + 1;

	tot_elem_num = Nx * Ny;
	tot_node_num = Nxl * Nyl;

	for (int i = 0; i < tot_node_num; i++) {
		Node n{};
		nodes.push_back(n);
	}

	for (int i = 0; i < tot_elem_num; i++) {
		Element e{};
		elements.push_back(e);
	}

	Eigen::VectorXd xn = Eigen::VectorXd::LinSpaced(Nxl, 0.0, L);
	Eigen::VectorXd yn = Eigen::VectorXd::LinSpaced(Nyl, 0.0, h);
	std::vector<std::vector<double>> node_position(tot_node_num, std::vector<double>(2, 0.0));

	for (size_t i = 0; i < Nxl; i++) {
		for (size_t j = 0; j < Nyl; j++) {
			node_position[static_cast<size_t>(i * Nyl + j)][0] = xn[i];
			node_position[static_cast<size_t>(i * Nyl + j)][1] = yn[j];
		}
	}

	for (int i = 0; i < tot_node_num; i++) {
		nodes[i].nid = i;
		nodes[i].gid[0] = (int)((nodes[i].nid + 1) * 2 - 1 - 1);
		nodes[i].gid[1] = (int)((nodes[i].nid + 1) * 2 - 1);
		nodes[i].position(0) = node_position[i][0];
		nodes[i].position(1) = node_position[i][1];
	}

	size_t i{};
	for (size_t e = 0; e < tot_elem_num; e++) {

		if ((e % Ny == 0) && (e != 0))
			i++;

		elements[e].eid = e;
		elements[e].n1 = &nodes[static_cast<size_t>(e + i)];
		elements[e].n2 = &nodes[static_cast<size_t>(e + Ny + 1 + i)];
		elements[e].n3 = &nodes[static_cast<size_t>(e + Ny + 2 + i)];
		elements[e].n4 = &nodes[static_cast<size_t>(e + 1 + i)];

	}		
}

Mesh::Mesh(const Mesh& msh) : tot_node_num{ msh.tot_node_num }, tot_elem_num{ msh.tot_elem_num },
							  elements{ msh.elements }, nodes{ msh.nodes } {}

Mesh::~Mesh() {}

