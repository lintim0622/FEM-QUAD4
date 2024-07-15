#pragma once

#include "calculate.h"

class Static_Solver {
public:

    Static_Solver() = delete;
    Static_Solver(Calculate& cal, std::vector<size_t>& bid_list);
    ~Static_Solver();

    // free dof matrix
    Eigen::MatrixXd M_ff;
    Eigen::MatrixXd K_ff;
    Eigen::MatrixXd F_f;

    std::vector<size_t> dof_id;
    size_t len_bid;

    size_t free_dof_num;
    size_t constraint_dof_num;
    size_t tot_dof_num;

    Eigen::MatrixXd global_U;
    // Eigen::MatrixXd U;
};