#include "calculate.h"
#include "shape_function.h"

#include <cmath>

Calculate::Calculate(Mesh& msh) : msh{ msh } {

    size_t dou_tot_node_num = this->msh.tot_node_num * 2;
    M = Eigen::MatrixXd::Zero(dou_tot_node_num, dou_tot_node_num);
    K = Eigen::MatrixXd::Zero(dou_tot_node_num, dou_tot_node_num);
    F = Eigen::MatrixXd::Zero(dou_tot_node_num, 1);

    pid = 0;
    p[0] = 0.0;
    p[1] = 0.0;

}

Calculate::~Calculate() {}

void Calculate::apply_force_info(size_t pid, double p[]) {

    this->pid = pid;
    this->p[0] = p[0];
    this->p[1] = p[1];

}

void Calculate::single_elem_matrix(Material& elastic) {

    // ''' strain-displacement matrix '''
    double& rho = elastic.rho();
    double& w = elastic.w();
    double& A = elastic.A();
    Eigen::Matrix<double, 3, 3>& D = elastic.D();
    
    Element& ieo = msh.elements[0];
    double& x1 = ieo.n1->position[0];
    double& y1 = ieo.n1->position[1];

    double& x2 = ieo.n2->position[0];
    double& y2 = ieo.n2->position[1];

    double& x3 = ieo.n3->position[0];
    double& y3 = ieo.n3->position[1];

    double& x4 = ieo.n4->position[0];
    double& y4 = ieo.n4->position[1];

    // w1 = w2 = 1.0 use 2 Gauss points
    Strain_Displacement eu{};
    Shape_Function sh{};

    double eta_arr[2] { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };
    double xi_arr[2] { -1.0 / std::sqrt(3.0), 1.0 / std::sqrt(3.0) };

    double x[4] { x1, x2, x3, x4 };
    double y[4] { y1, y2, y3, y4 };

    for (double &eta : eta_arr) {
        for (double &xi : xi_arr) {

            eu.local_dN(xi, eta);
            eu.Jacobian(x, y);
            eu.matrix();
            double& det_J = eu.det_J;
            Eigen::Matrix<double, 3, 8> B = eu.B;

            sh.local_N(xi, eta);
            sh.matrix();
            Eigen::Matrix<double, 2, 8>& N = sh.N;
            
            ieo.me += rho * A * N.transpose() * N * det_J;
            ieo.ke += (B.transpose() * D) * B * det_J * w;
        }
    }
}

void Calculate::apply_node_f_ext(const double& ti) {

    if (pid > msh.tot_node_num) {
        std::cerr << "reset external force id\n";
        return;
    }

    for (Node& node : msh.nodes) {
        node.f_ext.Zero();

        if (node.nid == pid && ti == 0.0) {
            node.f_ext(0) = this->p[0];
            node.f_ext(1) = this->p[1];
        }
    }
}

void Calculate::set_M_K_F() {

    // calculate K matrix and F vector
    Element& ieo = msh.elements[0];
    for (Element& ie : msh.elements) {
        double local_f[8]{ ie.n1->f_ext(0), ie.n1->f_ext(1), ie.n2->f_ext(0), ie.n2->f_ext(1),
                           ie.n3->f_ext(0), ie.n3->f_ext(1), ie.n4->f_ext(0), ie.n4->f_ext(1) };

        int gid[8]{ ie.n1->gid[0], ie.n1->gid[1], ie.n2->gid[0], ie.n2->gid[1],
                    ie.n3->gid[0], ie.n3->gid[1], ie.n4->gid[0], ie.n4->gid[1] };

        for (int i = 0; i < 8; i++) {
            F(gid[i], 0) = local_f[i];

            for (int j = 0; j < 8; j++) {
                M(gid[i], gid[j]) += ieo.me(i, j);
                K(gid[i], gid[j]) += ieo.ke(i, j);
            }
        }
    }

}