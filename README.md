# Sensor Fusion: Extended Kalman Filter Project

[//]: # (Image References)
[image1]: ./Final.png  "Final"
[image2]: ./KalmanFilterAlgo.png  "kalman"
[image3]: ./Missing_normalize.png  "Normalize"

The overall goal of this project is to utilize a kalman filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. The kalman filter is then run against a simulator to predict the motion of a car.  Lidar measurements are red circles and radar measurements are blue circles with an arrow pointing in the direction of the observed angle.  The path of the car is predetermined by the input data set, and the green triangles are the predicted path from the kalman filter.  Here is a snapshot of my final result:

![alt text][image1]

## [Rubric Points](https://review.udacity.com/#!/rubrics/748/view)
### Here I will consider the rubric points individually and describe how I addressed each point in my implementation.  

### Compiling

The main program can be built and run by doing the following from the project top directory.

1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./ExtendedKF

### Accuracy

The rubric specifies the maximum RMSE values for each measurement. My final RSME values vs. the limits were as follows:
| Measurement | Max RSME | Measured RSME  |
|:-----------:|:--------:|:--------------:|
|     px      |   0.11   |     0.0973     |
|     py      |   0.11   |     0.0855     |
|     vx      |   0.52   |     0.4513     |
|     vy      |   0.52   |     0.4399     |

Here's a [link to a video recording of my final result](./project_recording.mp4).

### Follows the Correct Algorithm
#### Sensor Fusion algorithm follows the general processing flow as taught in the preceding lessons

The Kalman filter I implemented generally follows the flow described in the lesson materials:

![alt text][image2]

#### Kalman Filter algorithm handles the first measurements appropriately

The 'FusionEKF' class initializes all necessary matrices and initializes a 'KalmanFilter' class to store these matrices and implement Update/Predict steps.  Input measurements enter the 'ProcessMeasurement()' function in the 'FusionEKF' class.  The measurements can be either Lidar or Radar data. On the first measurement, states are initialized based on the provided data.  For Radar data, the input polar coordinates are converted to cartesian coordinates using basic trigonometry and stored here.  

#### Kalman Filter algorithm first predicts then updates

On subsequent measurements, prediction is preformed using an updated covariance matrix and state transition matrix.  Then laser and/or radar matrices are setup, and states are updated with new measurement data. 

#### Kalman Filter can handle radar and lidar measurements

Lidar measurements use standard kalman filter equations to update states within the 'Update()' function of the 'KalmanFilter' class.  Standard measurement and measurement covariance matrices are used.  

Radar uses _extended_ kalman filter equations to update states within the 'UpdateEKF()' function.  For this, stored state data is  first converted to polar so that it can be compared against the new measurements (which are in polar coordinates).  Also, for an _extended_ kalman filter, the measurement matrix is a Jacobian matrix calculated in the 'CalculateJacobian()' function of the 'Tools' class.

### Code Efficiency
Algorithm does not need to sacrifice comprehension, stability, robustness or security for speed, while still maintaining good practice with respect to calculations (ex. avoided calling angle normalization in 'UpdateEKF()' unless angle was outside expected range, and avoided divide by zero in all cases).  Avoided code duplication where possible (ex. created 'UpdateHelper() function in 'KalmanFilter' class to avoid duplicate update/estimation code).

## Discussion
This project was less work intensive than previous ones as most of the framework was already setup in the sample code and in the previous lesson materials.  However, for Radar measurements in particular, there were a few details to work out that required me to understand the general concepts.  For example, conversion of the polar coordinates to cartesian in the initial measurement required going back to the basic trigonometry lesson.  

One thing that tripped me up was the normalization of the angle in y measurement.  This was mentioned in the "Tips and tricks" for the project, but I missed it during my first implementation.  Thus, my first attempt showed the following path when the turn angle switched from left to right at the middle of the test:

![alt text][image3]

I added some additional debugging to understand what was happening at this step.  I noticed my bearing measurement had switched from negative to positive (i.e. from turning left to turning right), which caused the difference between the last and current measurement to be larger than π.  As noticed in the "Tips and tricks", the Kalman filter expects small angle value differences between -π  and π .  Thus, I normalized the bearing portion of the "y" vector using atan2().  This fixed the issue at this point in the car's path, as seen in the final video.
