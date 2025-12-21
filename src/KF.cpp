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
