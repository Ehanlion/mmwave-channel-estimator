\# Background



This project was worked on for a Master's class of mine for UCLA. The project topic is about compressive estimation of mmWave channels. The goal of this project is to study the normalized mean square error (NMSE) of channel estimation as a function of SNR. Furthermore, we will study spectral efficiency as a function of the SNR. 



\# Core Concepts for Later



\## NMSE



NMSE = normalized mean square error. It is used to measure average difference between estimated and actual values. It is a good indicator of how accurate an estimation algorithm is. The NMSE is the normalized MSE, where MSE is the mean squared error. MSE is the average of squares of differences between prediction and actual. Normalizing this divides the MSE by variance of data set, making it scale-agnostic.



\## SNR



SNR compares the level of a desired signal against the level of noise. Basically, it measures the signal's clarity. Higher SNR means that the signal is clearer. Later, we plot NMSE vs. M and fix SNR to values. These values are -15, -5, 0 dB. 0dB here is the best case while -15 dB is more muddy. 



\## Spectral Efficiency



Spectral efficiency measures how efficiently a bandwidth is used to transmit data. Basically, it tells you how much data can be crammed into a given amount of bandwidth. Spectral efficiency (because I love cars) is like asking how many cars can I jam down a highway. The bandwidth is the number of lanes and the data is the cars travelling. My spectral efficiency is how many cars I can get moving here.



\# Reconstruction Algorithm (SW-OMP)



The main reconstruction task is to perform compressive recovery. This can be done with various algorithms but for this project, we are told to use Simultaneous Weighted Orthogonal Matching Pursuit, or SW-OMP. This is a greedy, iterative algorithm. It estimate one propagation path per iteration. After all estimations are done, it collects the gains of all paths estimated, refines them and their individual contributions, and then rectifies them against the original vector of measurements. It has a stop condition, which is basically just a number of iterations. In this project, we use M=80 or M=120 for our results. 



\# Deliverable Tasks



\## NMSE vs. SNR



Plot the NMSE vs. SNR over the range of -15 to 10 dB. Assume on-grid channel angles and generate two curves for the output, one for M=80 and M=120 (number of allowed iterations for the SW-OMP algorithm).



\## NMSE vs. Training Frames, SNR



Plot the NMSE vs. number of training frames. Vary M from 20 -> 100. Here we fix the SNR values to -10, -5, and 0 dB. So we will generate 3 plots, each with M varied over the range, one plot per fixed SNR value.



\## Spectral Efficiency vs. SNR



Plot the spectral efficiency vs. SNR over the range of -15 to 10 dB. For this task, assume on-grid channel angles and use a fixed number of training frames. Here we will fix M to just one value of 60. 



\# Project Layout



The following is the project layout. We will break the project down into subfiles for clarity. At the end we will create driver files for examples, to demonstrate and execute the 3 main deliverables for this project!



```txt

% source files

&nbsp; channel\_params.m          % struct of all constants (antennas, subcarriers, grids, etc.)

&nbsp; gen\_channel.m             % wideband mmWave channel (on-grid angles) per spec

&nbsp; build\_training.m          % training frames, RF/baseband training combining/precoding

&nbsp; vectorize\_measurements.m  % build y = A\*x per Sec. III-A (stack across frames)

&nbsp; swomp.m                   % Simultaneous Weighted OMP core (Alg. 1)

&nbsp; nmse.m                    

&nbsp; spectral\_efficiency.m     

&nbsp; seed\_rng.m                % reproducibility

% experiments

&nbsp; exp\_nmse\_vs\_snr.m

&nbsp; exp\_nmse\_vs\_M.m

&nbsp; exp\_rate\_vs\_snr.m

```



The project directory will also include the specs sheet for the project, an advisement paper about approaches and methods, as well as the final report. 

