# Assignment 11: Planning a Fast Rescue Boat Mission with an Extended Kalman Filter

### Group 6: Pavel Kuznetsov, Olga Tiupina, Sofya Konstantinova, Skoltech, 2025

## Code realization
### Kalman Filter class

EKF_11.hpp:

```cpp
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
```
EKF_11.cpp

```cpp
#include "EKF_11.hpp"

EKF::EKF(double x0x, double v0x,
         double x0y, double v0y,
         double varianceD, double varianceB,
         double varianceStateNoise,
         double timeStep, int N)
    : x0x(x0x), v0x(v0x), x0y(x0y), v0y(v0y),
      varianceD(varianceD),
      varianceB(varianceB),
      varianceStateNoise(varianceStateNoise),
      T(timeStep),
      N(N)
{
    initialization();
}

void EKF::initialization()
{

    // Prepare containers
    stateTrue.clear();
    statePredicted.clear();
    stateFiltered.clear();

    measurements.clear();

    Ppredicted.clear();
    Pfiltered.clear();

    stateTrue.reserve(N);
    statePredicted.reserve(N);
    stateFiltered.reserve(N);

    measurements.reserve(N);

    Pfiltered.reserve(N);
    Ppredicted.reserve(N);

    Eigen::Vector2d prev;
    prev << x0x, v0x;
    std::vector<Eigen::Vector2d> trajectoryX;
    generateTrueTrajectory(trajectoryX, prev);
    std::vector<Eigen::Vector2d> trajectoryY;
    prev << x0y, v0y;
    generateTrueTrajectory(trajectoryY, prev);
    mergeTrajectories(trajectoryX, trajectoryY);
    generateMeasurements();

    double x1m = measurements[0](0) * sin(measurements[0](1)); // t=2: D*sin(β)
    double y1m = measurements[0](0) * cos(measurements[0](1)); // t=2: D*cos(β)

    // initial state X0
    Eigen::Vector4d X0;

    X0 << x1m,
        0,
        y1m,
        0;
    stateFiltered.push_back(X0);
    statePredicted.push_back(X0);

    // initial P0,0
    Eigen::Matrix4d P0;
    P0 << std::pow(10, 10), 0.0, 0.0, 0.0,
        0.0, std::pow(10, 10), 0.0, 0.0,
        0.0, 0.0, std::pow(10, 10), 0.0,
        0.0, 0.0, 0.0, std::pow(10, 10);

    Pfiltered.push_back(P0);
    Ppredicted.push_back(P0);

    // R scalar for measurement
    R1 << varianceD, 0,
        0, varianceB;

    // Transition Phi and G
    Phi.setIdentity();
    Phi(0, 1) = T;
    Phi(2, 3) = T;
    G << (T * T) / 2.0, 0.0,
        T, 0.0,
        0.0, (T * T) / 2.0,
        0.0, T;

    // Q = G * G^T * variance_a
    Q = G * G.transpose() * varianceStateNoise;
}
void EKF::generateTrueTrajectory(std::vector<Eigen::Vector2d> &trajectory, Eigen::Vector2d prev)
{
    // Generate accelerations ~ N(meanAcc, varianceAcc)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> accDist(meanStateNoise, std::sqrt(varianceStateNoise));

    trajectory.push_back(prev);

    for (int k = 1; k < N; ++k)
    {
        double a = accDist(gen);
        double x_prev = trajectory[k - 1](0);
        double v_prev = trajectory[k - 1](1);

        double xk = x_prev + v_prev * T + 0.5 * a * T * T;
        double vk = v_prev + a * T;

        Eigen::Vector2d cur;
        cur << xk, vk;
        trajectory.push_back(cur);
    }
}

void EKF::mergeTrajectories(
    const std::vector<Eigen::Vector2d> &trajectoryX,
    const std::vector<Eigen::Vector2d> &trajectoryY)
{

    size_t size = trajectoryX.size();
    if (size != trajectoryY.size())
    {
        return;
    }

    stateTrue.clear();
    stateTrue.reserve(N);
    for (size_t i = 0; i < size; ++i)
    {
        stateTrue.emplace_back(
            trajectoryX[i](0), // x
            trajectoryX[i](1), // vx
            trajectoryY[i](0), // y
            trajectoryY[i](1)  // vy
        );
    }
}

void EKF::generateMeasurements()
{
    // Measurements noise ~ N(meanMeas, varianceMeas)
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> measDistD(0, std::sqrt(varianceD));
    std::normal_distribution<double> measDistB(0, std::sqrt(varianceB));

    measurements.clear();
    measurements.reserve(N);

    for (int k = 0; k < N; ++k)
    {

        double D_k = std::sqrt(std::pow(stateTrue[k](0), 2) + std::pow(stateTrue[k](2), 2));
        double b_k = std::atan2(stateTrue[k](0), stateTrue[k](2));
        D_k += measDistD(gen);
        b_k += measDistB(gen);
        Eigen::Vector2d measurement;
        measurement << D_k, b_k;
        measurements.push_back(measurement);
    }
}

// filter launch
void EKF::runFilter()
{
    // Already have X0|0 and P0|0 in stateFiltered[0], Pfiltered[0]
    for (int k = 1; k < N; ++k)
    {
        // Prediction
        Eigen::Vector4d Xprev = stateFiltered[k - 1];
        Eigen::Matrix4d Pprev = Pfiltered[k - 1];

        Eigen::Vector4d Xpred = Phi * Xprev;

        Eigen::Matrix4d Ppred = Phi * Pprev * Phi.transpose() + Q;

        statePredicted.push_back(Xpred);
        Ppredicted.push_back(Ppred);

        Eigen::Vector4d Xf;
        Eigen::Matrix4d Pf;
        Eigen::Matrix4d I = Eigen::Matrix4d::Identity();

        // radar measurement
        Eigen::Matrix<double, 2, 1> h_i;
        h_i << std::sqrt(Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2)), atan2(Xpred(0), Xpred(2));

        Eigen::Matrix<double, 2, 4> deltaH_i = Eigen::Matrix<double, 2, 4>::Zero();
        deltaH_i(0, 0) = Xpred(0) / std::sqrt(Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2));
        deltaH_i(0, 2) = Xpred(2) / std::sqrt(Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2));
        deltaH_i(1, 0) = Xpred(2) / (Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2));
        deltaH_i(1, 2) = (-1) * Xpred(0) / (Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2));

        Eigen::Vector2d z_k = measurements[k];
        Eigen::Matrix<double, 4, 2> K = Ppred * deltaH_i.transpose() * ((deltaH_i * Ppred * deltaH_i.transpose()) + R1).inverse();
        Xf = Xpred + K * (z_k - h_i);

        Pf = (I - K * deltaH_i) * Ppred;

        stateFiltered.push_back(Xf);
        Pfiltered.push_back(Pf);
    }



}

std::vector<Eigen::Vector4d> EKF::getFilteredState()
{
    return stateFiltered;
}

std::vector<Eigen::Vector4d> EKF::getPredictedState()
{
    return statePredicted;
}

std::vector<Eigen::Vector4d> EKF::getTrueState()
{
    return stateTrue;
}

std::vector<Eigen::Vector2d> EKF::getMeasuredState() {
    return measurements;
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
#include <iostream>

#include "Plotter.hpp"
#include "EKF_11.hpp"
#include "ErrorEstimater.hpp"



int main()
{
    // =========================================Trajectory generation and forward Kalman filter algorithm=========================================

    // Problem parameters per 1 part of assignment
    const int N = 500;
    const double T = 1.0;
    double varianceStateNoise = 0.3 * 0.3;
    double varianceD = 50 * 50;
    double varianceB = 0.004 * 0.004;
    const double x0x = 1000.0;
    const double x0y = 1000.0;
    const double v0x = 100.0;
    const double v0y = 100.0;

    // forward KF
    EKF ekf(x0x, v0x, x0y, v0y, varianceD, varianceB, varianceStateNoise, T, N);
    ekf.runFilter();

    std::vector<double> filtrationStepVec(N);
    for (int i = 0; i < N; i++) {
        filtrationStepVec[i] = i * T;
    }
    
    
       {
        std::vector<double> DTrue;
        std::vector<double> BTrue;
        auto statesTrue = ekf.getTrueState();
        for(const auto &stateTrue: statesTrue) {
            DTrue.push_back(std::sqrt(stateTrue(0)*stateTrue(0) + stateTrue(2)*stateTrue(2)));
            BTrue.push_back(std::atan2(stateTrue(0), stateTrue(2)));
        }

        auto statesFiltered = ekf.getFilteredState();
        
        std::vector<double> DFiltered;
        std::vector<double> BFiltered;
        for(const auto &stateFiltr: statesFiltered) {
            DFiltered.push_back(std::sqrt(stateFiltr(0)*stateFiltr(0) + stateFiltr(2)*stateFiltr(2)));
            BFiltered.push_back(std::atan2(stateFiltr(0), stateFiltr(2)));
        }
        
        auto statesPredicted = ekf.getPredictedState();
        std::vector<double> DPredicted;
        std::vector<double> BPredicted;
        for(const auto &statePred: statesPredicted) {
            DPredicted.push_back(std::sqrt(statePred(0)*statePred(0) + statePred(2)*statePred(2)));
            BPredicted.push_back(std::atan2(statePred(0), statePred(2)));
        }

        auto statesMeas = ekf.getMeasuredState();
        std::vector<double> DMeasured;
        std::vector<double> BMeasured;
        for(const auto &stateMeas: statesMeas) {
            DMeasured.push_back(stateMeas(0, 0));
            BMeasured.push_back(stateMeas(1, 0));
        }
        universalPlot({DTrue, DFiltered, DPredicted, DMeasured}, filtrationStepVec, {"True range, m", "Filtered estimate of range, m", "Predicted estimate of range, m", "Measurment of range, m"}, "Time, s", "Trajectories, m", "Fig. 1 Trajectories of range D");
        universalPlot({BTrue, BFiltered, BPredicted, BMeasured}, filtrationStepVec, {"True azimuth, rad", "Filtered estimate of azimuth, rad", "Predicted estimate of azimuth, rad", "Measurment of azimuth, rad"}, "Time, s", "Trajectories, rad", "Fig. 3 Trajectories of azimuth betta");
        
        const std::size_t crop = 220;
        auto cropBoth = [&](std::vector<double> &v) {
            if (v.size() <= 2 * crop) { v.clear(); return; }
            v.erase(v.begin(), v.begin() + crop);
            v.erase(v.end() - crop, v.end());
        };

        cropBoth(DTrue);
        cropBoth(DFiltered);
        cropBoth(DPredicted);
        cropBoth(DMeasured);

        cropBoth(BTrue);
        cropBoth(BFiltered);
        cropBoth(BPredicted);
        cropBoth(BMeasured);

        cropBoth(filtrationStepVec);

        universalPlot({DTrue, DFiltered, DPredicted, DMeasured}, filtrationStepVec, {"True range, m", "Filtered estimate of range, m", "Predicted estimate of range, m", "Measurment of range, m"}, "Time, s", "Trajectories, m", "Fig. 2 Cropped trajectories of range D");
        universalPlot({BTrue, BFiltered, BPredicted, BMeasured}, filtrationStepVec, {"True azimuth, rad", "Filtered estimate of azimuth, rad", "Predicted estimate of azimuth, rad", "Measurment of azimuth, rad"}, "Time, s", "Trajectories, rad", "Fig. 4 Cropped trajectories of azimuth betta");
    }

    // =========================================================Error===========================================================================

    std::vector<std::vector<double>> DFilteredVector;
    std::vector<std::vector<double>> BFilteredVector;
    std::vector<std::vector<double>> DPredictedVector;
    std::vector<std::vector<double>> BPredictedVector;
    std::vector<std::vector<double>> DTrueVector;
    std::vector<std::vector<double>> BTrueVector;

    for (int t = 0; t < 500; t++)
    {
        ekf.initialization();
        ekf.runFilter();

        auto statesFiltered = ekf.getFilteredState();
        std::vector<double> DFiltered;
        std::vector<double> BFiltered;
        for(const auto &stateFiltr: statesFiltered) {
            DFiltered.push_back(std::sqrt(stateFiltr(0)*stateFiltr(0) + stateFiltr(2)*stateFiltr(2)));
            BFiltered.push_back(std::atan2(stateFiltr(0), stateFiltr(2)));
        }
        DFilteredVector.push_back(DFiltered);
        BFilteredVector.push_back(BFiltered);

        auto statesPredicted = ekf.getPredictedState();
        std::vector<double> DPredicted;
        std::vector<double> BPredicted;
        for(const auto &statePred: statesPredicted) {
            DPredicted.push_back(std::sqrt(statePred(0)*statePred(0) + statePred(2)*statePred(2)));
            BPredicted.push_back(std::atan2(statePred(0), statePred(2)));
        }
        DPredictedVector.push_back(DPredicted);
        BPredictedVector.push_back(BPredicted);

        auto statesTrue = ekf.getTrueState();
        std::vector<double> DTrue;
        std::vector<double> BTrue;
        for(const auto &stateTrue: statesTrue) {
            DTrue.push_back(std::sqrt(stateTrue(0)*stateTrue(0) + stateTrue(2)*stateTrue(2)));
            BTrue.push_back(std::atan2(stateTrue(0), stateTrue(2)));
        }
        DTrueVector.push_back(DTrue);
        BTrueVector.push_back(BTrue);
    }

    std::vector<double> finalErrorXKFFilteredB;
    std::vector<double> finalErrorXKFFilteredD;
    std::vector<double> finalErrorXKFPredicteredB;
    std::vector<double> finalErrorXKFPredicteredD;


    ErrorEstimater::computeAllErrors(DTrueVector, DPredictedVector, finalErrorXKFPredicteredD);
    ErrorEstimater::computeAllErrors(DTrueVector, DFilteredVector, finalErrorXKFFilteredD);

    ErrorEstimater::computeAllErrors(BTrueVector, BPredictedVector, finalErrorXKFPredicteredB);
    ErrorEstimater::computeAllErrors(BTrueVector, BFilteredVector, finalErrorXKFFilteredB);

    std::cout << "finalErrorXKFPredicteredB: " << finalErrorXKFPredicteredB.size() << std::endl;
    std::cout << "finalErrorXKFFilteredB: " << finalErrorXKFFilteredB.size() << std::endl;
    
    std::vector<double> sigmaDNoise;
    sigmaDNoise.resize(N, std::sqrt(varianceD));

    std::vector<double> sigmaB1Noise;
    sigmaB1Noise.resize(N, std::sqrt(varianceB));

     // generate time steps for trajectories
    
    std::vector<double> errorStepVec(N - 3);
    for (int i = 3; i < N; i++) {
        errorStepVec[i - 3] = i * T;
    }
    universalPlot({finalErrorXKFPredicteredD, finalErrorXKFFilteredD, sigmaDNoise}, errorStepVec, {"True error KF for predicted estimate range", "True error KF for filtered estimate range", "Sigma noise D"}, "Time, s", "Error values, m", "Fig. 5 The true errors and std for the estimate of the rescue boat range");
    universalPlot({finalErrorXKFPredicteredB, finalErrorXKFFilteredB, sigmaB1Noise}, errorStepVec, {"True error KF for predicted estimate azimuth", "True error KF for filtered estimate azimuth", "Sigma noise B"}, "Time, s", "Error values, rad", "Fig. 6 The true errors and std for the estimate of the rescue boat azimuth");

    return 0;
}


```


# Results

![Alt text](<data_11/Fig. 1 Trajectories of range D.png>)

![Alt text](<data_11/Fig. 2 Cropped trajectories of range D.png>)

At each step of prediction and filtration, using an estimate of the predicted and filtered state vector X, the range D was calculated. The estimate of the predicted and filtered values is close to the true range value D throughout the entire observation time.

![Alt text](<data_11/Fig. 3 Trajectories of azimuth betta.png>)

![Alt text](<data_11/Fig. 4 Cropped trajectories of azimuth betta.png>)

At each step of prediction and filtration, using an estimate of the predicted and filtered state vector X, the azimuth beta was calculated. The estimate of the predicted and filtered values is close to the true azimuth beta value throughout the entire observation time.

![Alt text](<data_11/Fig. 5 The true errors and std for the estimate of the rescue boat range.png>)

Figure 5 shows that both the predicted and filtered range estimation errors quickly decrease  during the first 30-40 steps, and stabilize at a value around 15-17 m, which is significantly lower than the measurement noise level of 50 m.

![Alt text](<data_11/Fig. 6 The true errors and std for the estimate of the rescue boat azimuth.png>)

Figure 6 shows that both the predicted and filtered azimuth estimation errors quickly decrease during the first 30-40 steps to a value around 0.002 rad,  which is  2 times lower than the measurement noise level, and continues to decrease.

### Conclusion

During this assignment, we:

- To generate noisy radar measurements:

	- Generated the true trajectory of the boat as x and y coordinates, and calculated the true values of range D and azimuth beta

	- Then measurements of range D and azimuth beta were generated from the radar by adding noise to the true values of range B and azimuth beta

	- Next, the new x and y coordinates were calculated from the noisy measurements of the range D and azimuth beta from the radar

- Developed the Extended Kalman Filter (EKF) to estimate the position of a rescue boat moving along a nonlinear trajectory based on noisy radar measurements

- The developed EKF was used to obtain the predicted and filtered estimates of the range D and azimuth beta. These estimates closely followed the true trajectory.

- Analyzed the true estimation errors for both predicted and filtered estimation errors

	- For the range D estimation, the both errors quickly decreased during the first 30–40 steps and stabilized around 15–17 m, significantly below the measurement noise - 50 m.

	- For the azimuth beta estimation, the both errors quickly decrease during the first 30-40 steps to a value around 0.002 rad,  which is  2 times lower than the measurement noise level, and continues to decrease.

	- In both cases, the filtered estimation errors remained lower than the predicted estimation errors

Thus, the developed EKF provides an estimate of the range D and azimuth beta close to the real value based on noisy measurements from the radar throughout the entire observation time. The true errors of both estimates quickly decrease to values well below the measurement errors. This indicates that the EKF is correctly constructed to estimate the position in a nonlinear system.