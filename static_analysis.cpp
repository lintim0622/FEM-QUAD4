#include "static_analysis.h"

#include <algorithm>

Static_Solver::Static_Solver(Calculate& cal, std::vector<size_t>& bid_list)
    : free_dof_num {}, constraint_dof_num{}, len_bid{ bid_list.size() } {

    tot_dof_num = static_cast<size_t>(cal.msh.tot_node_num * 2);
    dof_id = std::vector<size_t>(tot_dof_num, 0);
    std::vector<size_t> nid_array(cal.msh.tot_node_num, 0);

    if (len_bid != 0) {

        std::sort(bid_list.begin(), bid_list.end());

        size_t i{};
        for (size_t& bid : bid_list) {
            nid_array[i] = bid;
            i++;
        }

        size_t m{ len_bid };
        for (size_t i{}; i < cal.msh.tot_node_num; i++) {
            size_t num{};

            for (size_t j{}; j < len_bid; j++) {
                if (i == nid_array[j])
                    num += (size_t)1;
            }

            if (num == 0) {
                nid_array[m] = i;
                m += static_cast<size_t>(1);
            }
        }

        m = (size_t)0;
        for (size_t& i : nid_array) {
            for (Node& node : cal.msh.nodes) {
                if (i == node.nid) {
                    dof_id[m] = static_cast<size_t>(node.gid[0]);
                    dof_id[static_cast<size_t>(m+1)] = static_cast<size_t>(node.gid[1]);
                    m += static_cast<size_t>(2);
                }
            }
        }

        free_dof_num = static_cast<size_t>(std::sqrt((tot_dof_num * tot_dof_num) - (len_bid * (size_t)2 * (tot_dof_num * (size_t)2 - len_bid * (size_t)2))));
        constraint_dof_num = static_cast<size_t>(tot_dof_num - free_dof_num);
        M_ff = Eigen::MatrixXd::Zero(free_dof_num, free_dof_num);
        K_ff = Eigen::MatrixXd::Zero(free_dof_num, free_dof_num);
        F_f = Eigen::MatrixXd::Zero(free_dof_num, 1);

        for (size_t i {}; i < free_dof_num; i++) {
            F_f(i, 0) = cal.F(dof_id[i + constraint_dof_num], 0);

            for (size_t j{}; j < free_dof_num; j++) {
                M_ff(i, j) = cal.M(dof_id[i + constraint_dof_num], dof_id[j + constraint_dof_num]);
                K_ff(i, j) = cal.K(dof_id[i + constraint_dof_num], dof_id[j + constraint_dof_num]);
            }
        }

        global_U = Eigen::MatrixXd::Zero(tot_dof_num, 1);
        Eigen::MatrixXd inv_K_ff = K_ff.inverse();
        Eigen::MatrixXd u_f = inv_K_ff * F_f;
        for (size_t i{}; i < free_dof_num; i++) {
            global_U(dof_id[i + constraint_dof_num], 0) = u_f(i, 0);
        }

        i = 0;
        for (Node& node : cal.msh.nodes) {
            node.displacement(0) = global_U(i, 0);
            node.displacement(1) = global_U(i+1, 0);
            i += 2;
        }
    }

    else {
        std::cerr << "No Boundary Condition\n";
        return;
    }
}

Static_Solver::~Static_Solver() {}