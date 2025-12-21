#include "EKF.hpp"

EKF::EKF(double x0x, double v0x,
         double x0y, double v0y,
         double varianceD, double varianceB,
         double varianceBAccur, double varianceStateNoise,
         double timeStep, int N)
    : x0x(x0x), v0x(v0x), x0y(x0y), v0y(v0y),
      varianceD(varianceD),
      varianceB(varianceB),
      varianceBAccur(varianceBAccur),
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

    
    double x1m = measurements[0](0) * sin(measurements[0](1));  // t=2: D*sin(β)
    double y1m = measurements[0](0) * cos(measurements[0](1));  // t=2: D*cos(β)
    double x3m = measurements[2](0) * sin(measurements[2](1));  // t=4: D*sin(β)  
    double y3m = measurements[2](0) * cos(measurements[2](1));  // t=4: D*cos(β)

    // initial state X0
    Eigen::Vector4d X0;
    
    X0 << x3m, 
          (x3m - x1m) / (2 * T),
          y3m,
          (y3m - y1m) / (2 * T);
    stateFiltered.push_back(X0);
    statePredicted.push_back(X0);

          // initial P0,0
    Eigen::Matrix4d P0;
    P0 << 10000, 0.0, 0.0, 0.0,
        0.0, 10000, 0.0, 0.0,
        0.0, 0.0, 10000, 0.0,
        0.0, 0.0, 0.0, 10000;
    Pfiltered.push_back(P0);
    Ppredicted.push_back(P0);

    // R scalar for measurement
    R1 << varianceD, 0,
        0, varianceB;
    R2 = varianceBAccur;

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
    std::normal_distribution<double> measDistBAccur(0, std::sqrt(varianceBAccur));

    measurements.clear();
    measurements.reserve(N);

    for (int k = 0; k < N; ++k)
    {

        if (k % 2 == 0)
        {
            double D_k = std::sqrt(std::pow(stateTrue[k](0), 2) + std::pow(stateTrue[k](2), 2));
            double b_k = std::atan2(stateTrue[k](0), stateTrue[k](2));
            D_k += measDistD(gen);
            b_k += measDistB(gen);
            Eigen::Vector2d measurement;
            measurement << D_k, b_k;
            measurements.push_back(measurement);
        }
        else
        {
            double D_k = std::numeric_limits<double>::quiet_NaN();
            double b_k = std::atan2(stateTrue[k](0), stateTrue[k](2));
            b_k += measDistBAccur(gen);
            Eigen::Vector2d measurement;
            measurement << D_k, b_k;
            measurements.push_back(measurement);
        }
    }
}

// filter launch
void EKF::runFilter()
{
    // Already have X0|0 and P0|0 in stateFiltered[0], Pfiltered[0]
    for (int k = 3; k < N; ++k) //since 4
    {
        // Prediction
        Eigen::Vector4d Xprev = stateFiltered[k - 3];
        Eigen::Matrix4d Pprev = Pfiltered[k - 3];

        Eigen::Vector4d Xpred = Phi * Xprev;
        Eigen::Matrix4d Ppred = Phi * Pprev * Phi.transpose() + Q;

        statePredicted.push_back(Xpred);
        Ppredicted.push_back(Ppred);

        Eigen::Vector4d Xf;
        Eigen::Matrix4d Pf;
        Eigen::Matrix4d I = Eigen::Matrix4d::Identity();

        // radar measurement
        if (k % 2 == 0)
        {
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
        }

        // NON-radar measurement
        else
        {
            Eigen::Matrix<double, 1, 1> h_i;
            h_i << atan2(Xpred(0), Xpred(2));

            Eigen::Matrix<double, 1, 4> deltaH_i = Eigen::Matrix<double, 1, 4>::Zero();
            deltaH_i(0, 0) = Xpred(2) / (Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2));
            deltaH_i(0, 2) = (-1) * Xpred(0) / (Xpred(0) * Xpred(0) + Xpred(2) * Xpred(2));

            double z_k = measurements[k](1);
            Eigen::Matrix<double, 1, 1> S = deltaH_i * Ppred * deltaH_i.transpose();
            S(0, 0) += R2;
            Eigen::Matrix<double, 4, 1> K = Ppred * deltaH_i.transpose() * S.inverse();
            Xf = Xpred + K * (z_k - h_i(0, 0));

            Pf = (I - K * deltaH_i) * Ppred;
        }

        stateFiltered.push_back(Xf);
        Pfiltered.push_back(Pf);
    }

}


std::vector<Eigen::Vector4d> EKF::getFilteredState() {
    return stateFiltered;
}

std::vector<Eigen::Vector4d> EKF::getPredictedState() {
    return statePredicted;
}

std::vector<Eigen::Vector4d> EKF::getTrueState() {
    return stateTrue;
}

std::vector<Eigen::Vector2d> EKF::getMeasuredState() {
    return measurements;
}
