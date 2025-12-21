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