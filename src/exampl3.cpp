#include <iostream>

#include "Plotter.hpp"
#include "EKF_11.hpp"
#include "ErrorEstimater.hpp"


std::vector<Eigen::Vector2d> readMeasurementsFromFile(const std::string &filename)
{
    std::vector<Eigen::Vector2d> measurements;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << "Error opening file: " << filename << std::endl;
        return measurements;
    }

    double first, second;
    while (file >> first >> second)
    {
        measurements.emplace_back(first, second);
    }

    file.close();
    return measurements;
}


int main()
{
    // =========================================Trajectory generation and forward Kalman filter algorithm=========================================

    std::vector<Eigen::Vector2d> z1 = readMeasurementsFromFile("/home/olia/Skoltech/DataProcessing/data_processing_25/final/data/z1.txt");
    std::vector<Eigen::Vector2d> z2 = readMeasurementsFromFile("/home/olia/Skoltech/DataProcessing/data_processing_25/final/data/z2.txt");
    std::vector<Eigen::Vector2d> z3 = readMeasurementsFromFile("/home/olia/Skoltech/DataProcessing/data_processing_25/final/data/z3.txt");

    // Problem parameters per 1 part of assignment
    const double T = 1.0;
    double varianceD = 200 * 200;
    double varianceB = 0.01 * 0.01;
    double varianceA = 0.01 * 0.01;


    double v0x = 0.0;
    double v0y = 0.0;


    double x0x1 = z1[0](0) * std::sin(z1[0](1));
    double x0y1 = z1[0](0) * std::cos(z1[0](1));

    double x0x2 = z2[0](0) * std::sin(z2[0](1));
    double x0y2 = z2[0](0) * std::cos(z2[0](1));

    double x0x3 = z3[0](0) * std::sin(z3[0](1));
    double x0y3 = z3[0](0) * std::cos(z3[0](1));

    // forward KF
    EKF ekf1(x0x1, v0x, x0y1, v0y, varianceD, varianceB, varianceA, T, z1);
    ekf1.runFilter();

    EKF ekf2(x0x2, v0x, x0y2, v0y, varianceD, varianceB, varianceA, T, z1);
    ekf2.runFilter();

    EKF ekf3(x0x3, v0x, x0y3, v0y, varianceD, varianceB, varianceA, T, z1);
    ekf3.runFilter();

    // smoothing KF
    
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

