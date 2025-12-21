#pragma once

#include <Eigen/Dense>
#include <vector>
#include <random>
#include <cmath>
#include <iostream>

/**
 * @class KalmanFilter
 * @brief Implementation of a Kalman Filter for 1D motion tracking
 *
 * This class implements a discrete Kalman Filter for tracking
 * position and velocity of an object with acceleration model.
 */
class KalmanFilter
{
public:
    KalmanFilter(double x0, double v0,
                 double varianceX0, double varianceV0,
                 double varianceAcc, double varianceMeas,
                 double timeStep, int N,
                 double meanAcc = 0.0, double meanMeas = 0.0)
        : x0(x0), v0(v0),
          varianceX0(varianceX0),
          varianceV0(varianceV0),
          varianceAcc(varianceAcc),
          varianceMeas(varianceMeas),
          meanAcc(meanAcc),
          meanMeas(meanMeas),
          T(timeStep),
          N(N)
    {
        initialization();
    }

    // filter launch
    void runFilter();

    void runSmoothing();

    // filter future prediction
    void predictFutureStates(int m);

    // initialization of all matrix
    void initialization();
    void initialization(double gapProbability);


    std::vector<double> getTrueTrajectoryX();
    std::vector<double> getTrueTrajectoryV();

    std::vector<double> getKalmanX();
    std::vector<double> getKalmanV();

    std::vector<double> getFutureX();
    std::vector<double> getFutureV();

    std::vector<double> getPredictX();

    std::vector<double> getMeasurments();

    std::vector<double> getSigmaX();
    std::vector<double> getSigmaV();

    std::vector<double> getSigmaPredX();
    std::vector<double> getSigmaFutureX();

    std::vector<double> getSigmaSmoothX();
    std::vector<double> getSigmaSmoothV();

    std::vector<double> getSmoothedX();
    std::vector<double> getSmoothedV();
    
    void futureError(int m);


private:
    void generateTrueTrajectoryAndMeasurements();
    void generateTrueTrajectoryAndGapMeasurements(double gapProbability);



    // model parameters / matrices
    double varianceAcc;
    double varianceMeas;
    double meanAcc;
    double meanMeas;
    double T;
    double x0;
    double v0;
    double varianceX0;
    double varianceV0;
    int N;

    Eigen::Matrix2d Phi;
    Eigen::Vector2d G;
    Eigen::Matrix2d Q;
    Eigen::Matrix<double, 1, 2> H;
    double R;

    std::vector<double> measurements; // z[k]

    std::vector<Eigen::Vector2d> stateTrue;      // X_true[k]
    std::vector<Eigen::Vector2d> stateFuture;    // X_{k + m - 1 |k}
    std::vector<Eigen::Vector2d> stateSmoothed;  // X_{k |N}
    std::vector<Eigen::Vector2d> statePredicted; // X_{k|k-1} stored at index k-1
    std::vector<Eigen::Vector2d> stateFiltered;  // X_{k|k}

    std::vector<Eigen::Matrix2d> Ppredicted; // P_{k|k-1}
    std::vector<Eigen::Matrix2d> Psmoothed;  // P_{k|N}
    std::vector<Eigen::Matrix2d> Pfiltered;  // P_{k|k}
    std::vector<Eigen::Matrix2d> Pfuture; 
    
    std::vector<Eigen::Matrix<double, 2, 1>> Kvec; // Kalman gains
    std::vector<Eigen::Matrix2d> Asmooth;          // Smoothed gain

    std::vector<double> sigmaX;       // sqrt(Pf(0,0))
    std::vector<double> sigmaV;       // sqrt(Pf(1,1))
    std::vector<double> sigmaPredX;   // sqrt(Pp(0,0))
    std::vector<double> sigmaPredV;   // sqrt(Pp(1,1))
    std::vector<double> sigmaSmoothX; // sqrt(Ps(0,0))
    std::vector<double> sigmaSmoothV; // sqrt(Ps(1,1))
    std::vector<double> sigmaFutureX; // sqrt(Pf(0,0))
    std::vector<double> sigmaFutureV; // sqrt(Pfut(1,1))
};
