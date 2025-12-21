# Assignment 12: Planning a Coastal Survey Aircraft Tracking Mission from Joint Assimilation of Heterogeneous Sensor data_12.

### Group 6: Pavel Kuznetsov, Olga Tiupina, Sofya Konstantinova, Skoltech, 2025

## Code realization
### Kalman Filter class

KF.hpp:

```cpp

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

```
KF.cpp

```cpp
#include "KF.hpp"

// filter launch
void KalmanFilter::runFilter()
{
    // Already have X0|0 and P0|0 in stateFiltered[0], Pfiltered[0]
    for (int k = 1; k < N; ++k)
    {
        // Prediction
        Eigen::Vector2d Xprev = stateFiltered[k - 1];
        Eigen::Matrix2d Pprev = Pfiltered[k - 1];

        Eigen::Vector2d Xpred = Phi * Xprev;
        Eigen::Matrix2d Ppred = Phi * Pprev * Phi.transpose() + Q;

        statePredicted.push_back(Xpred);
        Ppredicted.push_back(Ppred);
        sigmaPredX.push_back(std::sqrt(std::abs(Ppred(0, 0))));

        // Kalman gain: K = Ppred * H^T * inv(H*Ppred*H^T + R)
        Eigen::Matrix<double, 2, 1> K = Ppred * H.transpose() / ((H * Ppred * H.transpose())(0, 0) + R);
        Kvec.push_back(K);

        // Measurement at time k
        double z = measurements[k];
        Eigen::Vector2d Xf;
        Eigen::Matrix2d I;
        Eigen::Matrix2d Pf;

        if (std::isnan(z))
        {
            Xf = Xpred;
            Pf = Ppred;
        }
        else
        {
            // Filtration
            Xf = Xpred + K * (z - (H * Xpred)(0, 0));
            I = Eigen::Matrix2d::Identity();
            Pf = (I - K * H) * Ppred;
        }

        stateFiltered.push_back(Xf);
        Pfiltered.push_back(Pf);

        // Save sigma_x
        sigmaX.push_back(std::sqrt(std::abs(Pf(0, 0))));
        sigmaV.push_back(std::sqrt(std::abs(Pf(1, 1))));
    }
}

// filter future prediction
void KalmanFilter::predictFutureStates(int m)
{

    int n = static_cast<int>(stateFiltered.size());
    stateFuture.clear();
    // stateFuture.reserve(stateFiltered.size());
    for (int t = 0; t < n; t++)
    {
        stateFuture.push_back(Eigen::Vector2d::Zero());
    }

    Eigen::Matrix2d PhiFuture = Eigen::Matrix2d::Identity();
    for (int i = 0; i < m - 1; i++)
    {
        PhiFuture = PhiFuture * Phi;
    }

    for (int t = 0; t <= n - m; t++)
    {
        stateFuture[t + m - 1] = (PhiFuture * stateFiltered[t]);
    }
}

void KalmanFilter::futureError(int m)
{
    int n = static_cast<int>(stateFiltered.size());
    
    Pfuture.clear();
    sigmaFutureX.clear();
    sigmaFutureV.clear();

    Pfuture.resize(n, Eigen::Matrix2d::Zero());
    sigmaFutureX.resize(n, 0.0);
    sigmaFutureV.resize(n, 0.0);

    for (int i = 0; i < n - m; i++)
    {
        std::vector<Eigen::Matrix2d> PmSteps;
        PmSteps.push_back(Phi * Pfiltered[i] * Phi.transpose() + Q);

        for(int k = 1; k < m - 1; k++) {
            PmSteps.push_back(Phi * PmSteps[k - 1] * Phi.transpose() + Q);
        }
        
        Pfuture[i + m - 1] = PmSteps.back();

        sigmaFutureX[i + m - 1] = (std::sqrt(std::abs(PmSteps.back()(0, 0))));

        sigmaFutureV[i + m - 1] = (std::sqrt(std::abs(PmSteps.back()(1, 1))));

    }
}


// initialization of all matrix
void KalmanFilter::initialization()
{

    // Prepare containers
    stateTrue.clear();
    stateFuture.clear();
    statePredicted.clear();
    stateFiltered.clear();
    stateSmoothed.clear();

    measurements.clear();

    Psmoothed.clear();
    Ppredicted.clear();
    Pfiltered.clear();
    Pfuture.clear();

    Kvec.clear();
    Asmooth.clear();

    sigmaX.clear();
    sigmaPredX.clear();
    sigmaFutureX.clear();
    sigmaFutureV.clear();

    stateTrue.reserve(N);
    stateFuture.reserve(N);
    statePredicted.reserve(N);
    stateFiltered.reserve(N);
    stateSmoothed.reserve(N);

    measurements.reserve(N);

    Pfiltered.reserve(N);
    Ppredicted.reserve(N);
    Psmoothed.reserve(N);
    Pfuture.reserve(N);

    Kvec.reserve(N - 1);
    Asmooth.reserve(N);

    sigmaX.reserve(N);
    sigmaPredX.reserve(N);
    sigmaFutureX.reserve(N);
    sigmaFutureV.reserve(N);

    // initial state X0 = [x0; v0]
    Eigen::Vector2d X0;
    X0 << x0, v0;
    stateFiltered.push_back(X0);
    statePredicted.push_back(X0); // think about it

    // initial P0,0
    Eigen::Matrix2d P0;
    P0 << varianceX0, 0.0,
        0.0, varianceV0;
    Pfiltered.push_back(P0);
    Ppredicted.push_back(P0);

    // Observation matrix H = [1 0]
    H.setZero();
    H(0, 0) = 1.0;
    H(0, 1) = 0.0;

    // Transition Phi and G
    Phi.setIdentity();
    Phi(0, 1) = T;
    G << (T * T) / 2.0, T;

    // Q = G * G^T * variance_a
    Q = G * G.transpose() * varianceAcc;

    // R scalar for measurement
    R = varianceMeas;

    // set initial true state
    stateTrue.push_back(X0);

    // Generate the full true trajectory and measurements
    generateTrueTrajectoryAndMeasurements();
}

// initialization of all matrix
void KalmanFilter::initialization(double gapProbability)
{

    // Prepare containers
    stateTrue.clear();
    stateFuture.clear();
    statePredicted.clear();
    stateFiltered.clear();
    stateSmoothed.clear();

    measurements.clear();

    Psmoothed.clear();
    Ppredicted.clear();
    Pfiltered.clear();

    Kvec.clear();
    Asmooth.clear();

    sigmaX.clear();
    sigmaPredX.clear();

    stateTrue.reserve(N);
    stateFuture.reserve(N);
    statePredicted.reserve(N);
    stateFiltered.reserve(N);
    stateSmoothed.reserve(N);

    measurements.reserve(N);

    Pfiltered.reserve(N);
    Ppredicted.reserve(N);
    Psmoothed.reserve(N);

    Kvec.reserve(N - 1);
    Asmooth.reserve(N);

    sigmaX.reserve(N);
    sigmaPredX.reserve(N);

    // initial state X0 = [x0; v0]
    Eigen::Vector2d X0;
    X0 << x0, v0;
    stateFiltered.push_back(X0);
    statePredicted.push_back(X0); // think about it

    // initial P0,0
    Eigen::Matrix2d P0;
    P0 << varianceX0, 0.0,
        0.0, varianceV0;
    Pfiltered.push_back(P0);
    Ppredicted.push_back(P0);

    // Observation matrix H = [1 0]
    H.setZero();
    H(0, 0) = 1.0;
    H(0, 1) = 0.0;

    // Transition Phi and G
    Phi.setIdentity();
    Phi(0, 1) = T;
    G << (T * T) / 2.0, T;

    // Q = G * G^T * variance_a
    Q = G * G.transpose() * varianceAcc;

    // R scalar for measurement
    R = varianceMeas;

    // set initial true state
    stateTrue.push_back(X0);

    // Generate the full true trajectory and measurements
    generateTrueTrajectoryAndGapMeasurements(gapProbability);
}

void KalmanFilter::generateTrueTrajectoryAndMeasurements()
{
    // Generate accelerations ~ N(meanAcc, varianceAcc)
    // Measurements noise ~ N(meanMeas, varianceMeas)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> accDist(meanAcc, std::sqrt(varianceAcc));
    std::normal_distribution<double> measDist(meanMeas, std::sqrt(varianceMeas));

    for (int k = 1; k < N; ++k)
    {
        double a = accDist(gen);
        Eigen::Vector2d prev = stateTrue[k - 1];
        double x_prev = prev(0);
        double v_prev = prev(1);

        double xk = x_prev + v_prev * T + 0.5 * a * T * T;
        double vk = v_prev + a * T;

        Eigen::Vector2d cur;
        cur << xk, vk;
        stateTrue.push_back(cur);

        // measurement z_k = x_k + eta
        double z = xk + measDist(gen);
        measurements.push_back(z);
    }
    // measurement for t=0
    {
        double z0 = stateTrue[0](0) + measDist(gen);
        measurements.insert(measurements.begin(), z0);
    }
}

void KalmanFilter::generateTrueTrajectoryAndGapMeasurements(double gapProbability)
{
    // Generate accelerations ~ N(meanAcc, varianceAcc)
    // Measurements noise ~ N(meanMeas, varianceMeas)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> accDist(meanAcc, std::sqrt(varianceAcc));
    std::normal_distribution<double> measDist(meanMeas, std::sqrt(varianceMeas));

    // For generating gaps - uniform distribution for probability check
    std::uniform_real_distribution<double> gapDist(0.0, 1.0);
    // double gapProbability = 0.2; // P = 0.2

    for (int k = 1; k < N; ++k)
    {
        double a = accDist(gen);
        Eigen::Vector2d prev = stateTrue[k - 1];
        double x_prev = prev(0);
        double v_prev = prev(1);

        double xk = x_prev + v_prev * T + 0.5 * a * T * T;
        double vk = v_prev + a * T;

        Eigen::Vector2d cur;
        cur << xk, vk;
        stateTrue.push_back(cur);

        // Generate measurement with possible gap
        double xi = gapDist(gen);
        if (xi <= gapProbability)
        {
            measurements.push_back(std::numeric_limits<double>::quiet_NaN()); // Gap
        }
        else
        {
            double z = xk + measDist(gen);
            measurements.push_back(z);
        }
    }
    // measurement for t=0 (with possible gap)
    {
        double xi = gapDist(gen);
        if (xi <= gapProbability)
        {
            measurements.insert(measurements.begin(), std::numeric_limits<double>::quiet_NaN());
        }
        else
        {
            double z0 = stateTrue[0](0) + measDist(gen);
            measurements.insert(measurements.begin(), z0);
        }
    }
}

void KalmanFilter::runSmoothing()
{

    stateSmoothed.resize(N, Eigen::Vector2d::Zero());
    Asmooth.resize(N, Eigen::Matrix2d::Zero());
    Psmoothed.resize(N, Eigen::Matrix2d::Zero());

    sigmaSmoothX.resize(N - 1, 0.0);
    sigmaSmoothV.resize(N - 1, 0.0);

    // initialization for the last step
    stateSmoothed[N - 1] = stateFiltered[N - 1];
    Psmoothed[N - 1] = Pfiltered[N - 1];

    for (int t = N - 2; t >= 0; t--)
    {

        Asmooth[t] = Pfiltered[t] * Phi.transpose() * Ppredicted[t + 1].inverse();

        stateSmoothed[t] = stateFiltered[t] + Asmooth[t] * (stateSmoothed[t + 1] - Phi * stateFiltered[t]);
        Psmoothed[t] = Pfiltered[t] + Asmooth[t] * (Psmoothed[t + 1] - Ppredicted[t + 1]) * Asmooth[t].transpose();
        sigmaSmoothV[t] = std::sqrt(std::abs(Psmoothed[t](1, 1)));
        sigmaSmoothX[t] = std::sqrt(std::abs(Psmoothed[t](0, 0)));
    }
}

std::vector<double> KalmanFilter::getTrueTrajectoryX()
{
    std::vector<double> x(stateTrue.size());
    for (int t = 0; t < stateTrue.size(); t++)
    {
        x[t] = stateTrue[t](0);
    }
    return x;
}
std::vector<double> KalmanFilter::getTrueTrajectoryV()
{
    std::vector<double> v(stateTrue.size());
    for (int t = 0; t < stateTrue.size(); t++)
    {
        v[t] = stateTrue[t](1);
    }
    return v;
}

std::vector<double> KalmanFilter::getKalmanX()
{
    std::vector<double> x(stateFiltered.size());
    for (int t = 0; t < stateFiltered.size(); t++)
    {
        x[t] = stateFiltered[t](0);
    }
    return x;
}

std::vector<double> KalmanFilter::getKalmanV()
{
    std::vector<double> v(stateFiltered.size());
    for (int t = 0; t < stateFiltered.size(); t++)
    {
        v[t] = stateFiltered[t](1);
    }
    return v;
}

std::vector<double> KalmanFilter::getSmoothedX()
{
    std::vector<double> x(stateSmoothed.size());
    for (int t = 0; t < stateSmoothed.size(); t++)
    {
        x[t] = stateSmoothed[t](0);
    }
    return x;
}
std::vector<double> KalmanFilter::getSmoothedV()
{
    std::vector<double> v(stateSmoothed.size());
    for (int t = 0; t < stateSmoothed.size(); t++)
    {
        v[t] = stateSmoothed[t](1);
    }
    return v;
}

std::vector<double> KalmanFilter::getFutureX()
{
    std::vector<double> x(stateFuture.size());
    for (int t = 0; t < stateFuture.size(); t++)
    {
        x[t] = stateFuture[t](0);
    }
    return x;
}
std::vector<double> KalmanFilter::getFutureV()
{
    std::vector<double> v(stateFuture.size());
    for (int t = 0; t < stateFuture.size(); t++)
    {
        v[t] = stateFuture[t](1);
    }
    return v;
}
std::vector<double> KalmanFilter::getPredictX()
{
    std::vector<double> x(statePredicted.size());
    for (int t = 0; t < statePredicted.size(); t++)
    {
        x[t] = statePredicted[t](0);
    }

    return x;
}

std::vector<double> KalmanFilter::getMeasurments()
{
    return measurements;
}
std::vector<double> KalmanFilter::getSigmaX()
{
    return sigmaX;
}
std::vector<double> KalmanFilter::getSigmaV()
{
    return sigmaV;
}
std::vector<double> KalmanFilter::getSigmaPredX()
{
    return sigmaPredX;
}
std::vector<double> KalmanFilter::getSigmaFutureX()
{
    return sigmaFutureX;
}
std::vector<double> KalmanFilter::getSigmaSmoothX()
{
    return sigmaSmoothX;
}
std::vector<double> KalmanFilter::getSigmaSmoothV()
{
    return sigmaSmoothV;
}

```

### Visualization function realization

"Plotter.hpp"

```cpp
#pragma once

#include <matplot/matplot.h>
#include <string>

template<typename T>
void universalPlot(const std::vector<std::vector<double>> &axisesY, const std::vector<T> &axisX, std::vector<std::string> plotNames, std::string labelX, std::string labelY, std::string title) {
    static int num = 0;

    auto fig1 = matplot::figure();
    fig1->size(1200, 780);

    auto ax = matplot::gca();
    if(axisesY.size() > 1)
        matplot::hold(ax, true);
    
    int counter = 0;
    for(const auto axisY: axisesY) {
        matplot::plot(ax, axisX, axisY)->line_width(2).display_name(plotNames[counter]);
        counter++;
    }

    ax->xlabel(labelX);
    ax->ylabel(labelY);
    ax->title(title);
    ax->x_axis().label_font_size(18);
    ax->y_axis().label_font_size(18);
    ax->font_size(12);
    ax->title_font_size_multiplier(1.0);
    
    auto l = matplot::legend("show");
    l->font_size(16);
    l->location(matplot::legend::general_alignment::topright); 
    
    
    matplot::grid(true);
    matplot::save(title + ".png");
    matplot::show();
}
```

### Error finder 
ErrorEstimater.hpp

```cpp

#include <cmath>
#include <vector>

/**
 * @namespace ErrorEstimater
 * @brief Implementation of a two error types 
 * 
 * Calculate 
 * final error for each point for all measurments
 */
namespace ErrorEstimater {

    // put difference sqares to dynamic array (vector) 
    void errorI(const std::vector<double> &xTrue,
                const std::vector<double> &xEstimated,
                std::vector<double> &errorsVec) {
        if (xTrue.size() != xEstimated.size()) return;
        errorsVec.resize(xTrue.size());
        
        for (size_t t = 0; t < xTrue.size(); t++) {
            double diff = xTrue[t] - xEstimated[t];
            errorsVec[t] = diff * diff;
        }
    }
    // find errors
    void computeAllErrors(const std::vector<std::vector<double>> &xTrue,
                         const std::vector<std::vector<double>> &xEstimated,
                         std::vector<double> &finalError) {
        
        if (xTrue.size() != xEstimated.size() || xTrue.empty()) return;
        int KFlaunches = static_cast<int>(xTrue.size());
        int pointsPerLaunch = static_cast<int>(xTrue[0].size());
        
        std::vector<std::vector<double>> errorsVectors;
        errorsVectors.reserve(KFlaunches);
        
        for (int r = 0; r < KFlaunches; r++) {
            std::vector<double> errorsVector;
            errorI(xTrue[r], xEstimated[r], errorsVector);
            errorsVectors.push_back(std::move(errorsVector));
        }
        
        finalError.clear();
        finalError.reserve(pointsPerLaunch - 3);
        for (int i = 3; i < pointsPerLaunch; i++) {
            double meanErr = 0.0;
            for (int r = 0; r < KFlaunches; r++) {
                meanErr += errorsVectors[r][i];
            }
            meanErr /= KFlaunches;
            finalError.push_back(std::sqrt(meanErr));
        }
    }
}

```

### Main

```cpp

int main()
{
    // =========================================Trajectory generation and forward Kalman filter algorithm=========================================

    // Problem parameters per 1 part of assignment
    const int N = 200;
    const double T = 1.0;
    double varianceA = 0.2 * 0.2;
    const double varianceMeasNoise = 20.0 * 20.0;
    const double x0 = 5.0;
    const double v0 = 1.0;

    // initial uncertainty (suggested large values)
    double varX0 = 10000.0;
    double varV0 = 10000.0;

    // forward KF
    KalmanFilter kf(x0, v0, varX0, varV0, varianceA, varianceMeasNoise, T, N);
    kf.initialization(0.2);
    kf.runFilter();
    kf.predictFutureStates(7);

    // generate time steps for trajectories
    std::vector<int> filtrationStepVec(N - 1);
    for (auto &step : filtrationStepVec)
    {
        static int counter = 0;
        step = T * counter;
        counter++;
    }

    universalPlot({kf.getTrueTrajectoryX(), kf.getKalmanX(), kf.getMeasurments(), kf.getFutureX()}, filtrationStepVec, {"True trajectory", "Kalman filter trajectory", "Measurements", "Predicted future trajectory. m = 7"}, "Time, s", "Self-driving shuttle position, m", "Fig. 1 Estimate raft position");

    // =========================================================Error===========================================================================

    // P = 0.2

    std::vector<std::vector<double>> xTrue;
    std::vector<std::vector<double>> xKFFiltered;
    std::vector<std::vector<double>> xKFPredicted;
    std::vector<std::vector<double>> xKFFuture;

    for (int t = 0; t < 500; t++)
    {

        kf.initialization(0.2);
        kf.runFilter();
        kf.predictFutureStates(7);
        kf.futureError(7);

        xTrue.push_back(kf.getTrueTrajectoryX());
        xKFFiltered.push_back(kf.getKalmanX());
        xKFPredicted.push_back(kf.getPredictX());
        xKFFuture.push_back(kf.getFutureX());
    }

    std::vector<double> finalErrorXKFFutureP02;
    std::vector<double> finalErrorXKFFilteredP02;
    std::vector<double> finalErrorXKFPredicteredP02;

    ErrorEstimater::computeAllErrors(xTrue, xKFFuture, finalErrorXKFFutureP02);
    ErrorEstimater::computeAllErrors(xTrue, xKFFiltered, finalErrorXKFFilteredP02);
    ErrorEstimater::computeAllErrors(xTrue, xKFPredicted, finalErrorXKFPredicteredP02);

    std::vector<double> calculErrorXKFFilteredP02 = kf.getSigmaX();
    std::vector<double> calculErrorXKFFutureP02 = kf.getSigmaFutureX();
    std::vector<double> calculErrorXKFPredicteredP02 = kf.getSigmaPredX();

    universalPlot({finalErrorXKFFutureP02, kf.getSigmaX()}, filtrationStepVec, {"True error KF for filtered estimate, X", "Calculated error KF for filtered estimate, X"}, "Time, s", "Error values, m", "Fig. 2 The true errors and the calculated errors for the filtered estimates of the coordinate");
    universalPlot({finalErrorXKFFilteredP02, kf.getSigmaFutureX()}, filtrationStepVec, {"True error KF for predicted future estimates for m = 7, X", "Calculated error KF for predicted future estimates for m = 7, X"}, "Time, s", "Error values, m", "Fig. 3 The true errors and the calculated errors for the predicted future estimates for m = 7 of the coordinate");
    universalPlot({finalErrorXKFPredicteredP02, kf.getSigmaPredX()}, filtrationStepVec, {"True error KF for predicted estimates, X", "Calculated error KF for predicted estimates, X"}, "Time, s", "Error values, m", "Fig. 4 The true errors and the calculated errors for the predicted estimates of the coordinate");

    // P = 0.3

    xTrue.clear();
    xKFFiltered.clear();
    xKFPredicted.clear();
    xKFFuture.clear();

    for (int t = 0; t < 500; t++)
    {

        kf.initialization(0.3);
        kf.runFilter();
        kf.predictFutureStates(7);
        kf.futureError(7);

        xTrue.push_back(kf.getTrueTrajectoryX());
        xKFFiltered.push_back(kf.getKalmanX());
        xKFPredicted.push_back(kf.getPredictX());
        xKFFuture.push_back(kf.getFutureX());
    }

    std::vector<double> finalErrorXKFFutureP03;
    std::vector<double> finalErrorXKFFilteredP03;
    std::vector<double> finalErrorXKFPredicteredP03;

    ErrorEstimater::computeAllErrors(xTrue, xKFFuture, finalErrorXKFFutureP03);
    ErrorEstimater::computeAllErrors(xTrue, xKFFiltered, finalErrorXKFFilteredP03);
    ErrorEstimater::computeAllErrors(xTrue, xKFPredicted, finalErrorXKFPredicteredP03);

    std::vector<double> calculErrorXKFFilteredP03 = kf.getSigmaX();
    std::vector<double> calculErrorXKFFutureP03 = kf.getSigmaFutureX();
    std::vector<double> calculErrorXKFPredicteredP03 = kf.getSigmaPredX();


    // P = 0.5

    xTrue.clear();
    xKFFiltered.clear();
    xKFPredicted.clear();
    xKFFuture.clear();

    for (int t = 0; t < 500; t++)
    {

        kf.initialization(0.5);
        kf.runFilter();
        kf.predictFutureStates(7);
        kf.futureError(7);

        xTrue.push_back(kf.getTrueTrajectoryX());
        xKFFiltered.push_back(kf.getKalmanX());
        xKFPredicted.push_back(kf.getPredictX());
        xKFFuture.push_back(kf.getFutureX());
    }

    std::vector<double> finalErrorXKFFutureP05;
    std::vector<double> finalErrorXKFFilteredP05;
    std::vector<double> finalErrorXKFPredicteredP05;

    ErrorEstimater::computeAllErrors(xTrue, xKFFuture, finalErrorXKFFutureP05);
    ErrorEstimater::computeAllErrors(xTrue, xKFFiltered, finalErrorXKFFilteredP05);
    ErrorEstimater::computeAllErrors(xTrue, xKFPredicted, finalErrorXKFPredicteredP05);

    std::vector<double> calculErrorXKFFilteredP05 = kf.getSigmaX();
    std::vector<double> calculErrorXKFFutureP05 = kf.getSigmaFutureX();
    std::vector<double> calculErrorXKFPredicteredP05 = kf.getSigmaPredX();

    // P = 0.7

    xTrue.clear();
    xKFFiltered.clear();
    xKFPredicted.clear();
    xKFFuture.clear();

    for (int t = 0; t < 500; t++)
    {

        kf.initialization(0.7);
        kf.runFilter();
        kf.predictFutureStates(7);
        kf.futureError(7);

        xTrue.push_back(kf.getTrueTrajectoryX());
        xKFFiltered.push_back(kf.getKalmanX());
        xKFPredicted.push_back(kf.getPredictX());
        xKFFuture.push_back(kf.getFutureX());
    }

    std::vector<double> finalErrorXKFFutureP07;
    std::vector<double> finalErrorXKFFilteredP07;
    std::vector<double> finalErrorXKFPredicteredP07;

    ErrorEstimater::computeAllErrors(xTrue, xKFFuture, finalErrorXKFFutureP07);
    ErrorEstimater::computeAllErrors(xTrue, xKFFiltered, finalErrorXKFFilteredP07);
    ErrorEstimater::computeAllErrors(xTrue, xKFPredicted, finalErrorXKFPredicteredP07);

    std::vector<double> calculErrorXKFFilteredP07 = kf.getSigmaX();
    std::vector<double> calculErrorXKFFutureP07 = kf.getSigmaFutureX();
    std::vector<double> calculErrorXKFPredicteredP07 = kf.getSigmaPredX();


    universalPlot({finalErrorXKFFilteredP02, finalErrorXKFFilteredP03, finalErrorXKFFilteredP05, finalErrorXKFFilteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "True error KF for filtered estimate x, m", "Fig. 5 The final errors for the filtered estimates of the coordinate");
    universalPlot({finalErrorXKFPredicteredP02, finalErrorXKFPredicteredP03, finalErrorXKFPredicteredP05, finalErrorXKFPredicteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "True error KF for filtered estimate x, m", "Fig. 6 The final errors for the predicted estimates of the coordinate");
    universalPlot({finalErrorXKFFutureP02, finalErrorXKFFutureP03, finalErrorXKFFutureP05, finalErrorXKFFutureP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "True error KF for future estimate (m = 7) x, m", "Fig. 7 The final errors for the future predicted estimates (m = 7) of the coordinate");


    universalPlot({calculErrorXKFFilteredP02, calculErrorXKFFilteredP03, calculErrorXKFFilteredP05, calculErrorXKFFilteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for filtered estimate x, m", "Fig. 8 The calculated errors for the filtered estimates of the coordinate");
    universalPlot({calculErrorXKFPredicteredP02, calculErrorXKFPredicteredP03, calculErrorXKFPredicteredP05, calculErrorXKFPredicteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for filtered estimate x, m", "Fig. 9 The calculated errors for the predicted estimates of the coordinate");
    universalPlot({calculErrorXKFFutureP02, calculErrorXKFFutureP03, calculErrorXKFFutureP05, calculErrorXKFFutureP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for future estimate (m = 7) x, m", "Fig. 10 The calculated errors for the future predicted estimates (m = 7) of the coordinate");

    calculErrorXKFFilteredP02.erase(calculErrorXKFFilteredP02.begin(), calculErrorXKFFilteredP02.begin() + 20);
    calculErrorXKFFilteredP03.erase(calculErrorXKFFilteredP03.begin(), calculErrorXKFFilteredP03.begin() + 20);
    calculErrorXKFFilteredP05.erase(calculErrorXKFFilteredP05.begin(), calculErrorXKFFilteredP05.begin() + 20);
    calculErrorXKFFilteredP07.erase(calculErrorXKFFilteredP07.begin(), calculErrorXKFFilteredP07.begin() + 20);
    filtrationStepVec.erase(filtrationStepVec.begin(), filtrationStepVec.begin() + 10);
    universalPlot({calculErrorXKFFilteredP02, calculErrorXKFFilteredP03, calculErrorXKFFilteredP05, calculErrorXKFFilteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for filtered estimate x, m", "Fig. 11 The calculated errors for the filtered estimates of the coordinate");
    
    calculErrorXKFPredicteredP02.erase(calculErrorXKFPredicteredP02.begin(), calculErrorXKFPredicteredP02.begin() + 20);
    calculErrorXKFPredicteredP03.erase(calculErrorXKFPredicteredP03.begin(), calculErrorXKFPredicteredP03.begin() + 20);
    calculErrorXKFPredicteredP05.erase(calculErrorXKFPredicteredP05.begin(), calculErrorXKFPredicteredP05.begin() + 20);
    calculErrorXKFPredicteredP07.erase(calculErrorXKFPredicteredP07.begin(), calculErrorXKFPredicteredP07.begin() + 20);
    universalPlot({calculErrorXKFPredicteredP02, calculErrorXKFPredicteredP03, calculErrorXKFPredicteredP05, calculErrorXKFPredicteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for filtered estimate x, m", "Fig. 12 The calculated errors for the predicted estimates of the coordinate");
    
    calculErrorXKFFutureP02.erase(calculErrorXKFFutureP02.begin(), calculErrorXKFFutureP02.begin() + 20);
    calculErrorXKFFutureP03.erase(calculErrorXKFFutureP03.begin(), calculErrorXKFFutureP03.begin() + 20);
    calculErrorXKFFutureP05.erase(calculErrorXKFFutureP05.begin(), calculErrorXKFFutureP05.begin() + 20);
    calculErrorXKFFutureP07.erase(calculErrorXKFFutureP07.begin(), calculErrorXKFFutureP07.begin() + 20);
    universalPlot({calculErrorXKFFutureP02, calculErrorXKFFutureP03, calculErrorXKFFutureP05, calculErrorXKFFutureP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for future estimate (m = 7) x, m", "Fig. 13 The calculated errors for the future predicted estimates (m = 7) of the coordinate");

    return 0;
}

```


# Results

![Alt text](<data_12/Fig. 1 Trajectories of range D.png>)

![Alt text](<data_12/Fig. 2 Cropped trajectories of range D.png>)

At each step of prediction and filtration, using an estimate of the predicted and ltered state vector X, the range D was calculated. The estimate of the predicted and ltered values is close to the true range value D throughout the entire observation time.

![Alt text](<data_12/Fig. 3 Trajectories of azimuth betta.png>)

![Alt text](<data_12/Fig. 4 Cropped trajectories of azimuth betta.png>)

At each step of prediction and filtration, using an estimate of the predicted and filtered state vector X, the azimuth beta was calculated. 

Since the measurements of the azimuth beta alternately come from two different observers, when calculating the estimate of the filtered value of the azimuth beta:
- the measurement vector Z
- the observation function h(X)
- the measurement noise covariance matrix R

took different forms depending on the observer.

The estimate of the predicted and filtered values is close to the true azimuth beta value throughout the entire observation time.

![Alt text](<data_12/Fig. 5 The true errors and std for the estimate of the range.png>)

Figure 5 shows that both the predicted estimation and filtered range estimation errors quickly decrease during the first 30-40 steps, and stabilize at a value around 27-33 m, which approximately is two times lower than the measurement noise level of 50 m.

The true error of the prediction estimation and the filtered estimation range values D have maximums and minimums at different steps:
- **The maximum** true error of the **predicted** estimates is reached for **odd steps**, because these estimates were predicted on the **even steps**, when we had not information about range D
- **The minimum** true error of the **predicted** estimates is reached for **even steps**, because these estimates were predicted on the **odd steps**, when we had information about range D
- **The maximum** true error of the **filtered** estimates is reached on the **even steps**, because these estimates were filtered on the **even steps**, when we had not information about range D
- **The minimum** true error of the **filtered** estimates is reached on the **odd steps**, because these estimates were filtered on the **odd steps**, when we had information obout range D

![Alt text](<data_12/Fig. 6 The true errors and std for the estimate of the azimuth.png>)

Figure 6 shows that both the predicted estimation and filtered azimuth estimation errors quickly decreased  during the first 30-40 steps to a value of about 0.001 rad, which is 4 times lower than the noise level of the first observer and approximately equal to the noise level of the second observer, then the error continues to decrease.

The true error of the prediction estimation and the filtered azimuth beta estimation have maximums and minimums at different steps:
- **The maximum** true error of the **predicted** estimates is achieved for **even steps**, because these estimates were predicted for **odd steps**, when we received measurements of the azimuth beta from the first observer (less accurate)
- **The minimum** true error of the **predicted** estimates is achieved for **odd steps**, because these estimates were predicted for **even steps**, when we received measurements of the azimuth beta from the second observer (more accurate)
- **The maximum** true error of the **filtered** estimates is achieved on the **odd steps**, because these estimates were predicted for **even steps**, when we received measurements of the azimuth beta from the first observer (less accurate)
- **The minimum** true error of the **filtered** estimates is achieved on the **even steps**, because these estimates were predicted for **odd steps**, when we received measurements of the azimuth beta from the second observer (more accurate)

### Conclusion

During this assignment, we:
- Generated radar measurements with interference:
	- We generated the true trajectory of the boat in the form of x and y coordinates and calculated the true values of range D and azimuth beta
	- The range D and beta azimuth measurements were then obtained using radar by adding noise to the true range B and beta azimuth values
	- Moreover, at odd steps we receive information about the range D and azimuth beta from the observer with lower accuracy, and at even steps we receive information only about azimuth beta, but with better accuracy
	- To simulate receiving data from different observers when switching to a new filtering step, there was a switch between different ones: the measurement vector Z, the observation function h(X) and the measurement error covariance matrix R
	- The new x and y coordinates were then calculated based on noisy measurements of the D range and beta azimuth using radar

- An advanced Kalman filter (EKF) has been developed to estimate the location of a rescue boat moving along a nonlinear trajectory based on radar measurements with interference from two observers with different measurement accuracy
- The developed EKF was used to obtain predicted and filtered estimates of the D range and beta azimuth. These estimates exactly corresponded to the true trajectory

- Analyzed the true estimation errors for both predicted and ltered estimation errors:
	- For the range D estimation, the both errors quickly decreased during the first 30–40 steps and stabilized at a value around 27–33 m, which approximately is two times lower than the measurement noise level of 50 m.
	- For the azimuth beta estimation, the both errors quickly decrease during the first 30-40 steps to a value of about 0.001 rad, which is 4 times lower than the noise level of the first observer and approximately equal to the noise level of the second observer, then the error continues to decrease
	- In both cases, the ltered estimation errors remained lower than the predicted estimation errors

- The true error of the prediction estimation and the filtered estimation range values D and azimuth beta have maximums and minimums at different steps:
	- For range D:
		- **The maximum** true error of the **predicted** estimates is reached for **odd steps**
		- **The minimum** true error of the **predicted** estimates is reached for **even steps**
		- **The maximum** true error of the **filtered** estimates is reached on the **even steps**
		- **The minimum** true error of the **filtered** estimates is reached on the **odd steps**

	- For azimuth beta:
		- **The maximum** true error of the **predicted** estimates is achieved for **even steps**
		- **The minimum** true error of the **predicted** estimates is achieved for **odd steps**
		- **The maximum** true error of the **filtered** estimates is achieved on the **odd steps**
		- **The minimum** true error of the **filtered** estimates is achieved on the **even steps**

Thus, the developed EKF provides an estimate of the range D and azimuth beta close to the real value based on noisy measurements from two observers with different accuracy of the measurements throughout the entire observation time. The true errors of both estimates quickly decrease to values well below the measurement errors. This indicates that the EKF is correctly constructed to estimate the position in a nonlinear system.