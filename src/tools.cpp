#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {
    residue_ = VectorXd(4);
    residue_ << 0, 0, 0, 0;
    last_estimation_size_ = 0;
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

    VectorXd rmse(4);
    rmse << 0,0,0,0;

    for(int i=last_estimation_size_; i < estimations.size(); ++i){
        VectorXd err = estimations[i] - ground_truth[i];
        err = err.array() * err.array();
        residue_ += err;
    }
    last_estimation_size_ = estimations.size();

    rmse = residue_ / estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;

}