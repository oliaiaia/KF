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