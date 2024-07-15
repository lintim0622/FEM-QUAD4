#pragma once

#include "mesh.h"

class Calculate {
public:

    Calculate() = delete;
    Calculate(Mesh& msh);
    ~Calculate();

    void apply_force_info(size_t pid, double p[]);
    void single_elem_matrix(Material& elastic);
    void apply_node_f_ext(const double& ti);
    void set_M_K_F();

    // model
    Mesh& msh;

    // global mass matrix
    Eigen::MatrixXd M;

    // global stiffness matrix
    Eigen::MatrixXd K;

    // global force vector
    Eigen::MatrixXd F;

    // force information
    size_t pid;
    double p[2];
};