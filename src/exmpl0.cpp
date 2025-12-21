#include <iostream>

#include "Plotter.hpp"
#include "KF.hpp"
#include "ErrorEstimater.hpp"

std::vector<double> readNumbersFromFile(const std::string &filename)
{
    std::vector<double> data;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        std::cerr << filename << std::endl;
        return data;
    }

    double value;
    while (file >> value)
    {
        data.push_back(value);
    }

    file.close();
    return data;
}

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

    universalPlot({finalErrorXKFFilteredP02, kf.getSigmaX()}, filtrationStepVec, {"True error KF for filtered estimate, X", "Calculated error KF for filtered estimate, X"}, "Time, s", "Error values, m", "Fig. 2 The true errors and the calculated errors for the filtered estimates of the coordinate");
    universalPlot({finalErrorXKFFutureP02, kf.getSigmaFutureX()}, filtrationStepVec, {"True error KF for predicted future estimates for m = 7, X", "Calculated error KF for predicted future estimates for m = 7, X"}, "Time, s", "Error values, m", "Fig. 3 The true errors and the calculated errors for the predicted future estimates for m = 7 of the coordinate");
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


    calculErrorXKFFilteredP02.erase(calculErrorXKFFilteredP02.begin(), calculErrorXKFFilteredP02.begin() + 10);
    calculErrorXKFFilteredP03.erase(calculErrorXKFFilteredP03.begin(), calculErrorXKFFilteredP03.begin() + 10);
    calculErrorXKFFilteredP05.erase(calculErrorXKFFilteredP05.begin(), calculErrorXKFFilteredP05.begin() + 10);
    calculErrorXKFFilteredP07.erase(calculErrorXKFFilteredP07.begin(), calculErrorXKFFilteredP07.begin() + 10);
    filtrationStepVec.erase(filtrationStepVec.begin(), filtrationStepVec.begin() + 10);
    universalPlot({calculErrorXKFFilteredP02, calculErrorXKFFilteredP03, calculErrorXKFFilteredP05, calculErrorXKFFilteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for filtered estimate x, m", "Fig. 8 The calculated errors for the filtered estimates of the coordinate");
    
    calculErrorXKFPredicteredP02.erase(calculErrorXKFPredicteredP02.begin(), calculErrorXKFPredicteredP02.begin() + 10);
    calculErrorXKFPredicteredP03.erase(calculErrorXKFPredicteredP03.begin(), calculErrorXKFPredicteredP03.begin() + 10);
    calculErrorXKFPredicteredP05.erase(calculErrorXKFPredicteredP05.begin(), calculErrorXKFPredicteredP05.begin() + 10);
    calculErrorXKFPredicteredP07.erase(calculErrorXKFPredicteredP07.begin(), calculErrorXKFPredicteredP07.begin() + 10);
    universalPlot({calculErrorXKFPredicteredP02, calculErrorXKFPredicteredP03, calculErrorXKFPredicteredP05, calculErrorXKFPredicteredP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for filtered estimate x, m", "Fig. 9 The calculated errors for the predicted estimates of the coordinate");
    
    calculErrorXKFFutureP02.erase(calculErrorXKFFutureP02.begin(), calculErrorXKFFutureP02.begin() + 10);
    calculErrorXKFFutureP03.erase(calculErrorXKFFutureP03.begin(), calculErrorXKFFutureP03.begin() + 10);
    calculErrorXKFFutureP05.erase(calculErrorXKFFutureP05.begin(), calculErrorXKFFutureP05.begin() + 10);
    calculErrorXKFFutureP07.erase(calculErrorXKFFutureP07.begin(), calculErrorXKFFutureP07.begin() + 10);
    universalPlot({calculErrorXKFFutureP02, calculErrorXKFFutureP03, calculErrorXKFFutureP05, calculErrorXKFFutureP07}, filtrationStepVec, {"P = 0.2", "P = 0.3", "P = 0.5", "P = 0.7"}, "Time, s", "Calculated error KF for future estimate (m = 7) x, m", "Fig. 10 The calculated errors for the future predicted estimates (m = 7) of the coordinate");

    return 0;
}
