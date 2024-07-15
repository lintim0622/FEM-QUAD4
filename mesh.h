#ifndef MESH_H
#define MESH_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <cstdlib>

class Material {
public:
	Material();
    Material(double rho, double E, double v);
	~Material();

	void inputData();
    void cross_section(double h, double w);
	void outputFile();

    // get function
    double& rho();
    double& w();
    double& A();
    Eigen::Matrix<double, 3, 3>& D();
    

private:
	double _rho;
	double E;
	double v;

    // cross section
    double h;
    double _w;
    double _A;
    double I;

    // plane stress
    double E1;
    double E2;
    double G;
    Eigen::Matrix<double, 3, 3> _D;
};

class Node {
public:

    // function
    Node();
    ~Node();

    // member
    int nid;
    int gid[2];
    Eigen::Matrix<double, 2, 1> position;
    Eigen::Matrix<double, 2, 1> displacement;
    Eigen::Matrix<double, 3, 1> stress;
    Eigen::Matrix<double, 3, 1> strain;
    Eigen::Matrix<double, 2, 1> f_ext;
};

class Element {
public:

    // function
    Element();
    ~Element();

    // member
    int eid;
    Node* n1;
    Node* n2;
    Node* n3;
    Node* n4;
    Eigen::Matrix<double, 8, 8> me;
    Eigen::Matrix<double, 8, 8> ke;
};

class Mesh {
public:

    // function
    Mesh() = delete;
    Mesh(double L, double h, size_t Nx, size_t Ny);
    Mesh(const Mesh& msh);
    ~Mesh();

    // member
    size_t tot_node_num;
    size_t tot_elem_num;
    std::vector<Element> elements;
    std::vector<Node> nodes;
};

#endif // !MESH_H
