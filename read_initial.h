#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>


bool readTxtFile(const std::string& filename, std::map<std::string, std::vector<double>>& variables);

void apply_parameter(std::map<std::string, std::vector<double>>& variables,
    double& L, double& h, double& w, size_t& Nx, size_t& Ny,
    double& rho, double& E, double& v, int& pid, double P[],
    std::vector<size_t>& bid_list);