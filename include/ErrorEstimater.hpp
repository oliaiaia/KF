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

    // basic checks
    if (xTrue.empty() || xEstimated.empty()) return;

    // allow different numbers of KF launches (use the minimal available)
    int KFlaunches = static_cast<int>(std::min(xTrue.size(), xEstimated.size()));

        // determine minimal points-per-launch across all launches to avoid OOB
        size_t minPoints = std::numeric_limits<size_t>::max();
        for (int r = 0; r < KFlaunches; r++) {
            // safe access only within KFlaunches
            size_t sTrue = xTrue[r].size();
            size_t sEst = xEstimated[r].size();
            minPoints = std::min(minPoints, std::min(sTrue, sEst));
        }
        if (minPoints == std::numeric_limits<size_t>::max() || minPoints <= 3) return;

        // compute squared errors per launch up to minPoints
        std::vector<std::vector<double>> errorsVectors;
        errorsVectors.reserve(KFlaunches);

        for (int r = 0; r < KFlaunches; r++) {
            std::vector<double> errorsVector(minPoints);
            for (size_t t = 0; t < minPoints; t++) {
                double diff = xTrue[r][t] - xEstimated[r][t];
                errorsVector[t] = diff * diff;
            }
            errorsVectors.push_back(std::move(errorsVector));
        }

        // aggregate mean sqrt error starting from index 3
        finalError.clear();
        finalError.reserve(minPoints - 3);
        for (size_t i = 3; i < minPoints; i++) {
            double meanErr = 0.0;
            for (int r = 0; r < KFlaunches; r++) {
                meanErr += errorsVectors[r][i];
            }
            meanErr /= static_cast<double>(KFlaunches);
            finalError.push_back(std::sqrt(meanErr));
        }
    }
}

