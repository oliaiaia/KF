#include <iostream>

#include "Plotter.hpp"
#include "EKF.hpp"
#include "ErrorEstimater.hpp"



int main()
{
    // =========================================Trajectory generation and forward Kalman filter algorithm=========================================

    // Problem parameters per 1 part of assignment
    const int N = 500;
    const double T = 2.0;
    double varianceStateNoise = 0.3 * 0.3;
    double varianceD = 50 * 50;
    double varianceB = 0.004 * 0.004;
    double varianceBAccur = 0.001 * 0.001;
    const double x0x = 1000.0;
    const double x0y = 1000.0;
    const double v0x = 100.0;
    const double v0y = 100.0;

    // forward KF
    EKF ekf(x0x, v0x, x0y, v0y, varianceD, varianceB, varianceBAccur, varianceStateNoise, T, N);
    ekf.initialization();
    ekf.runFilter();

    // generate time steps for trajectories
    
    std::vector<double> filtrationStepVec(N - 2);
    // fill times for the filtered/predicted estimate indices (0..N-3) -> times (2*T .. (N-1)*T)
    for (int i = 0; i < N - 2; i++) {
        filtrationStepVec[i] = (i + 2) * T;
    }
    
    {
        std::vector<double> DTrue;
        std::vector<double> BTrue;
        auto statesTrue = ekf.getTrueState();
        for(const auto &stateTrue: statesTrue) {
            DTrue.push_back(std::sqrt(stateTrue(0)*stateTrue(0) + stateTrue(2)*stateTrue(2)));
            BTrue.push_back(std::atan2(stateTrue(0), stateTrue(2)));
        }
        DTrue.erase(DTrue.begin(), DTrue.begin() + 2);
        BTrue.erase(BTrue.begin(), BTrue.begin() + 2);

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
            if(!std::isnan(stateMeas(0))) {
                DMeasured.push_back(stateMeas(0, 0));
            }
            BMeasured.push_back(stateMeas(1, 0));
        }
        DMeasured.erase(DMeasured.begin(), DMeasured.begin() + 1);
        BMeasured.erase(BMeasured.begin(), BMeasured.begin() + 2);

        std::vector<double> measStepVec;
        // fill times for the filtered/predicted estimate indices (0..N-3) -> times (2*T .. (N-1)*T)
        for (int i = 0; i < N - 2; i+=2) {
            measStepVec.push_back((i + 2) * T);
        }

        std::vector<std::vector<double>> DSeries = {DTrue, DFiltered, DPredicted, DMeasured};
        std::vector<std::vector<double>> stepSeries = {filtrationStepVec, filtrationStepVec, filtrationStepVec, measStepVec};
        std::vector<std::string> DLabels = {"True range, m", "Filtered estimate of range, m", "Predicted estimate of range, m", "Measurment of range, m"};
        universalPlotVec(DSeries, stepSeries, DLabels, "Time, s", "Trajectories, m", "Fig. 1 Trajectories of range D");
        
        std::vector<std::vector<double>> BSeries = {BTrue, BFiltered, BPredicted, BMeasured};
        std::vector<std::string> BLabels = {"True azimuth, rad", "Filtered estimate of azimuth, rad", "Predicted estimate of azimuth, rad", "Measurment of azimuth, rad"};
        universalPlot(BSeries, filtrationStepVec, BLabels, "Time, s", "Trajectories, rad", "Fig. 3 Trajectories of azimuth betta");
        
        const std::size_t crop = 240;
        auto cropBoth = [&](std::vector<double> &v) {
            if (v.size() <= 2 * crop) { v.clear(); return; }
            v.erase(v.begin(), v.begin() + crop);
            v.erase(v.end() - crop, v.end());
        };

        cropBoth(DTrue);
        cropBoth(DFiltered);
        cropBoth(DPredicted);

        cropBoth(BTrue);
        cropBoth(BFiltered);
        cropBoth(BPredicted);
        cropBoth(BMeasured);

        cropBoth(filtrationStepVec);

        measStepVec.erase(measStepVec.begin(), measStepVec.begin() + 120);
        measStepVec.erase(measStepVec.end() - 120, measStepVec.end());

        DMeasured.erase(DMeasured.begin(), DMeasured.begin() + 120);
        DMeasured.erase(DMeasured.end() - 120, DMeasured.end());
        
        DSeries = {DTrue, DFiltered, DPredicted, DMeasured};
        stepSeries = {filtrationStepVec, filtrationStepVec, filtrationStepVec, measStepVec};
        DLabels = {"True range, m", "Filtered estimate of range, m", "Predicted estimate of range, m", "Measurment of range, m"};
        universalPlotVec(DSeries, stepSeries, DLabels, "Time, s", "Trajectories, m", "Fig. 2 Cropped trajectories of range D");
        
        
        universalPlot({BTrue, BFiltered, BPredicted, BMeasured}, filtrationStepVec, {"True azimuth, rad", "Filtered estimate of azimuth, rad", "Predicted estimate of azimuth, rad", "Measurment of azimuth, rad"}, "Time, s", "Trajectories, rad", "Fig. 4 Cropped trajectories of azimuth betta");
    }

    // // =========================================================Error===========================================================================

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
        //because start from k = 3
        DTrue.erase(DTrue.begin(), DTrue.begin() + 2);
        BTrue.erase(BTrue.begin(), BTrue.begin() + 2);

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
    
    std::vector<double> sigmaDNoise;
    sigmaDNoise.resize(N - 2, std::sqrt(varianceD));

    std::vector<double> sigmaB1Noise;
    sigmaB1Noise.resize(N - 2, std::sqrt(varianceB));
    std::vector<double> sigmaB2Noise;
    sigmaB2Noise.resize(N - 2, std::sqrt(varianceBAccur));

    std::vector<double> errorStepVec(N - 2);
    // fill times for the filtered/predicted estimate indices (0..N-3) -> times (2*T .. (N-1)*T)
    for (int i = 0; i < N - 2; i++) {
        errorStepVec[i] = (i + 2) * T;
    }
    universalPlot({finalErrorXKFPredicteredD, finalErrorXKFFilteredD, sigmaDNoise}, errorStepVec, {"True error KF for predicted estimate range", "True error KF for filtered estimate range", "Sigma noise D"}, "Time, s", "Error values, m", "Fig. 5 The true errors and std for the estimate of the range");
    universalPlot({finalErrorXKFPredicteredB, finalErrorXKFFilteredB, sigmaB1Noise, sigmaB2Noise}, errorStepVec, {"True error KF for predicted estimate azimuth", "True error KF for filtered estimate azimuth", "Sigma noise B", "Sigma noise more accurate B"}, "Time, s", "Error values, rad", "Fig. 6 The true errors and std for the estimate of the azimuth");


    return 0;
}
