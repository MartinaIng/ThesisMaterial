function [mean_error,var_error] = error_analysis(error)
% function that computes the mean and the variance of the errors
% 
% INPUT
% error: vector with the error for each time
% 
% OUTPUTS
% mean_error: mean of the error
% var_error: variance of the error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Computation of some statistical measures for the error time series
mean_error = mean(error);
var_error = var(error);
end