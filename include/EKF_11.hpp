#pragma once

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

class EKF
{
public:
    EKF(double x0x, double v0x,
             double x0y, double v0y,
                 double varianceD, double varianceB,
                double varianceStateNoise,
                 double timeStep, int N);
    void runFilter();
    void initialization();
    std::vector<Eigen::Vector4d> getFilteredState();
    std::vector<Eigen::Vector4d> getPredictedState();
    std::vector<Eigen::Vector4d> getTrueState();
    std::vector<Eigen::Vector2d> getMeasuredState();
    
    private:

    void mergeTrajectories(const std::vector<Eigen::Vector2d>& trajectoryX, const std::vector<Eigen::Vector2d>& trajectoryY);
    void generateTrueTrajectory(std::vector<Eigen::Vector2d> &trajectory, Eigen::Vector2d prev);
    void generateMeasurements();


    std::vector<Eigen::Vector4d> stateTrue;      // X_true[k]
    std::vector<Eigen::Vector4d> statePredicted; // X_{k|k-1} stored at index k-1
    std::vector<Eigen::Vector4d> stateFiltered;  // X_{k|k}

    std::vector<Eigen::Vector2d> measurements; // |D, betta|
    
    std::vector<Eigen::Matrix4d> Ppredicted; // P_{k|k-1}
    std::vector<Eigen::Matrix4d> Pfiltered;  // P_{k|k}


    Eigen::Matrix4d Phi;
    Eigen::Matrix<double, 4, 2> G;
    Eigen::Matrix4d Q;
    Eigen::Matrix2d R1;

    double varianceD, varianceB;
    double meanStateNoise = 0, varianceStateNoise;
    double x0x, v0x, x0y, v0y;

    int T, N;
};