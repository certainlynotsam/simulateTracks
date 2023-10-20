% simulateTracks.m
% written by Sam Daly 

clear all; close all; clc; 

addpath('lib')

%% GENERATE TRACKS

% Parameters
motion_type = ['brownian']; % 'brownian' / 'confined' / 'mixed'
space_units = 'µm'; % dimension units
time_units = 's'; % exposure time units

N_particles = 1; % number of particles
N_time_steps = 4000; % track lengths
N_dim = 3; % number of dimensions

% Area size, just used to disperse particles in 2D (no impact on analysis.)
SIZE = 10; % µm

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
Diff  = 0.1; % µm^2/s

% Time step between each acquisition
dT = 0.02; % s,

tracks = Types.generateTracks(motion_type, Diff,dT,N_particles,N_time_steps,N_dim,SIZE); % Generate tracks

Plotting.plotTracks(tracks, N_particles,space_units,1) % Plot tracks

% %% Measuring Diffusion with MLE
% 
% % Initialise a structure to hold the key bits of data of importance. 
% D = struct('TrackNo', zeros(N_particles,1), 'TrackDataPoints', zeros(N_particles,1), ...
%     'DxyzTrans', zeros(N_particles,1), 'DxyzTransError', zeros(N_particles,1), 'DxTrans', ...
%     zeros(N_particles,1), 'DxTransError', zeros(N_particles,1), 'DyTrans', zeros(N_particles,1), ...
%     'DyTransError',  zeros(N_particles,1), 'DzTrans', zeros(N_particles,1), 'DzTransError', zeros(N_particles,1));
% 
% % Run analysis of diffusion using MLE over all tracks
% for i = 1:max(N_particles)
% 
%     data = tracks{i};
% 
%     MLEF = MLEs; 
%         
% %% Calculate translational diffusion parameters using MLE        
% % Get trajectory coordinates
% t = data(1:end,1);
% TimeStep = diff(t(1:2));
% X = data(1:end,2);
% Y = data(1:end,3);
% Z = data(1:end,4);
% 
% [TrackDataPoints,RRTrackX,DtrueX, ... 
%     SigmaX,SigmaXError,RRTrackY, ...
%     DtrueY,SigmaY,SigmaYError, ...
%     RRTrackZ,DtrueZ,SigmaZ, ...
%     SigmaZError,MinXIndex,MinYIndex, ...
%     MinZIndex,DxTrans,DxTransError, ... 
%     DyTrans,DyTransError,DzTrans, ...
%     DzTransError,DxyzTrans,DxyzTransError,DcorrX,DcorrY,DcorrZ,N] = ...
%     MLEF.transMLE(X,Y,Z,TimeStep,length(X)-1);   
% 
% if DxyzTrans > 0 && DxTrans > 0 && DyTrans > 0 && DzTrans > 0
%     
%     D.TrackNo(i,1) = i;
%     D.TrackDataPoints(i,1) = TrackDataPoints;
%     D.DxyzTrans(i,1) = DxyzTrans*1e12;
%     D.DxyzTransError(i,1) = DxyzTransError*1e12;
%     D.DxTrans(i,1) = DxTrans*1e12;
%     D.DxTransError(i,1) = DxTransError*1e12;
%     D.DyTrans(i,1) = DyTrans*1e12;
%     D.DyTransError(i,1) = DyTransError*1e12;
%     D.DzTrans(i,1) = DzTrans*1e12;
%     D.DzTransError(i,1) = DzTransError*1e12;
% 
%     else 
% 
%     D.TrackNo(i,1) = NaN;
%     D.TrackDataPoints(i,1) = NaN;
%     D.DxyzTrans(i,1) = NaN;
%     D.DxyzTransError(i,1) = NaN;
%     D.DxTrans(i,1) = NaN;
%     D.DxTransError(i,1) = NaN;
%     D.DyTrans(i,1) = NaN;
%     D.DyTransError(i,1) = NaN;
%     D.DzTrans(i,1) = NaN;
%     D.DzTransError(i,1) = NaN;
% 
%     end
%  
%     progress = i/N_particles*100;
%     TF = isinteger(int8(progress));
%     if TF == 1
%     fprintf('MLE progress: %.1f%%\n',progress)
%     else end 
% 
% end
% 
% 
%     D.TrackNo = rmmissing(D.TrackNo);
%     D.TrackDataPoints = rmmissing(D.TrackDataPoints);
%     D.DxyzTrans = rmmissing(D.DxyzTrans);
%     D.DxyzTransError = rmmissing(D.DxyzTransError);
%     D.DxTrans = rmmissing(D.DxTrans);
%     D.DxTransError = rmmissing(D.DxTransError);
%     D.DyTrans = rmmissing(D.DyTrans);
%     D.DyTransError = rmmissing(D.DyTransError);
%     D.DzTrans = rmmissing(D.DzTrans);
%     D.DzTransError = rmmissing(D.DzTransError);
% 
% 
% % median diffusion coefficient in xyz 
% fprintf('D_{MLE}: %.5f +/- %.5f\n',median(D.DxyzTrans),median(D.DxyzTransError))
% 
% % histogram of diffusion coefficient in xyz
%     figure(2)
%     histogram(D.DxyzTrans,'BinWidth',0.01);
%     hold on
%     xline(Diff)
%     hold off
%     ylabel('Frequency');
%     xlabel('{\it D}_{xyz} (µm s^{-1})');

%% Measuring diffusion with MSD

MSD_fit_percent = 0.2;
loglog_fit = 0.6;

ma = msdanalyzer(N_dim, space_units, time_units);
ma = ma.addAll(tracks);

% Calculate MSD over all tracks
ma = ma.computeMSD;
figure
hmsd = ma.plotMeanMSD(gca, true);

[fo, gof] = ma.fitMeanMSD( MSD_fit_percent );
plot(fo)
legend off
ma.labelPlotMSD

% Fit individual MSD curves and calculate then mean using first 25%
ma = ma.fitMSD;

good_enough_fit = ma.lfit.r2fit > 0.8;
Dmean = mean( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;
Dstd  =  std( ma.lfit.a(good_enough_fit) ) / 2 / ma.n_dim;

fprintf('Estimation of the diffusion coefficient from linear fit of the 25%% MSD curves:\n')
fprintf('D = %.3g ± %.3g (mean ± std, N = %d)\n', Dmean, Dstd, sum(good_enough_fit));

% Determine the mean value of a log-log fit to determine confinement
ma = ma.fitLogLogMSD( loglog_fit );
ma.loglogfit;

LOGmean = mean(ma.loglogfit.alpha(good_enough_fit));
LOGstd = std(ma.loglogfit.alpha(good_enough_fit));

fprintf('Estimation of alpha from linear fit of %.0f %% of curve:\n', loglog_fit*100)
fprintf('a = %.3g ± %.3g (mean ± std, N = %d)\n', LOGmean, LOGstd, sum(good_enough_fit));


