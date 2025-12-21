# Results

## Measurement generation with gap = 0.2

The trajectories were ploted:

![alt text](<data/Fig. 1 Estimate raft position.png>)

The measurement data contains deliberate gaps in accordance with the problem statement. Predictive estimates for future time steps demonstrate higher variance and increased susceptibility to outliers relative to the filtered estimates.

By plot of errors for gap with probability equals to 0.2:


For filtered estimate:

![alt text](<data/Fig. 2 The true errors and the calculated errors for the filtered estimates of the coordinate.png>)


For predicted estimate:

![alt text](<data/Fig. 3 The true errors and the calculated errors for the predicted future estimates for m = 7 of the coordinate.png>)


For future predicted estimate (m = 7):

![alt text](<data/Fig. 4 The true errors and the calculated errors for the predicted estimates of the coordinate.png>)

For each estimation type, we observe close alignment between the true errors and calculated theoretical errors. This correspondence serves as evidence of filter optimality and proper parameter tuning.

## Comparing errors with different gap

### True error

The true error for filtered estimate:

![alt text](<data/Fig. 5 The final errors for the filtered estimates of the coordinate.png>)


The true error for predicted estimate:

![alt text](<data/Fig. 6 The final errors for the predicted estimates of the coordinate.png>)


The true error for predicted future estimate:

![alt text](<data/Fig. 7 The final errors for the future predicted estimates (m = 7) of the coordinate.png>)

The analysis reveals a clear dependence between the initial covariance P and estimation accuracy. For P = 0.2, the steady-state error reaches the lowest value, demonstrating superior estimation performance. As P increases to 0.7, the estimation error grows significantly, indicating degraded filter performance with overestimated initial uncertainty.

### Calculated error

This dependence is not so obvious, so we cut the beginning:

The calculated for filtered estimate:

![alt text](<data/Fig. 8 The calculated errors for the filtered estimates of the coordinate.png>)


The calculated for predicted estimate:

![alt text](<data/Fig. 9 The calculated errors for the predicted estimates of the coordinate.png>)


The calculated for predicted future estimate:

![alt text](<data/Fig. 10 The calculated errors for the future predicted estimates (m = 7) of the coordinate.png>)


All of them showed the same relationship. Calculated errors increase with increasing gap probability.

This report presents a comprehensive analysis of Kalman filter performance for urban shuttle tracking in challenging GNSS environments characterized by frequent signal dropouts. The study demonstrates successful implementation of a robust tracking system capable of maintaining accurate position estimates and generating reliable multi-step predictions despite significant measurement gaps.

1. Trajectory Generation System
- True Trajectory was generated using kinematic equations with normally distributed random acceleration
- Gap Simulation was implemented probabilistic measurement dropouts using uniform random variable ξ with varying probabilities P = [0.2, 0.3, 0.5, 0.7]

2. Kalman Filter Implementation
Developed an Kalman filter with handling for measurement gaps.

Main results:

1. Error Analysis with Gap Probability P = 0.2
The analysis demonstrated filter optimality: true errors and calculated errors are close across all estimation types.

2. Impact of Gap Probability on Estimation Accuracy
All errors increas with gap probability increasing.

This study successfully demonstrates that Kalman filter performance in urban shuttle tracking exhibits a clear inverse relationship with measurement availability. Both true errors and calculated errors increase with higher gap probabilities, confirming the expectation that reduced observational data degrades estimation accuracy.