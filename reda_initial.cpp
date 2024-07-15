#include "read_initial.h"

bool readTxtFile(const std::string& filename, std::map<std::string, std::vector<double>>& variables) {

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file " << "\"" << filename << "\"" << std::endl;
        std::cout << "This cpp file path¡G" << __FILE__ << std::endl;
        return false;
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#')
            continue;

        std::istringstream iss(line);
        std::string variableName;
        if (!(iss >> variableName))
            continue;

        if (variableName == "P") {
            std::vector<double> values;
            double value;
            while (iss >> value) {
                values.push_back(value);
            }
            if (values.size() >= 2) {
                variables["P"] = values;
            }
        }

        else if (variableName == "bid_list") {
            std::vector<double> values;
            double value;
            while (iss >> value) {
                values.push_back(value);
            }

            if (!values.empty()) {
                variables["bid_list"] = values;
            }
        }

        else {
            double value;
            if (iss >> value)
                variables[variableName].push_back(value);
        }
    }

    file.close();
    return true;
}

void apply_parameter(std::map<std::string, std::vector<double>>& variables,
    double& L, double& h, double& w, size_t& Nx, size_t& Ny,
    double& rho, double& E, double& v, int& pid, double P[],
    std::vector<size_t>& bid_list) {

    L = variables["L"][0];
    h = variables["h"][0];
    w = variables["w"][0];
    std::cout << "L: " << L << std::endl;
    std::cout << "h: " << h << std::endl;
    std::cout << "w: " << w << std::endl;
    std::cout << std::endl;

    Nx = variables["Nx"][0];
    Ny = variables["Ny"][0];
    std::cout << "Nx: " << Nx << std::endl;
    std::cout << "Ny: " << Ny << std::endl;
    std::cout << std::endl;

    rho = variables["rho"][0];
    E = variables["E"][0];
    v = variables["v"][0];
    std::cout << "rho: " << rho << std::endl;
    std::cout << "E: " << E << std::endl;
    std::cout << "v: " << v << std::endl;
    std::cout << std::endl;

    pid = variables["pid"][0];
    if (variables.count("pid") > 0) {
        std::cout << "pid: " << pid << std::endl;
    }

    else {
        std::cerr << "Error: pid not found in the file." << std::endl;
    }

    if (variables.count("P") > 0) {
        std::vector<double> tm_P = variables["P"];
        P[0] = tm_P[0];
        P[1] = tm_P[1];
        std::cout << "P:";
        for (double force : tm_P) {
            std::cout << " " << force;
        }
        std::cout << std::endl;
    }

    else {
        std::cerr << "Error: P not found in the file." << std::endl;
    }

    if (variables.count("bid_list") > 0) {
        std::vector<double> tmp_bid = variables["bid_list"];
        bid_list.push_back(static_cast<size_t>(tmp_bid[0]));
        bid_list.push_back(static_cast<size_t>(tmp_bid[1]));
        std::cout << "bid_list:";
        for (double id : bid_list) {
            std::cout << " " << id;
        }
        std::cout << std::endl;
    }

    else {
        std::cerr << "Error: bid_list not found in the file." << std::endl;
    }
}