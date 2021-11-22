function [alpha,gamma,eps] = lsfit(deg,y,rho)
% function that computes the parameters of the probability distribution function
% (using least square approximation) and the mean fitting error
% 
% INPUTS
% deg: vector with the degree for each node
% y: vector with the probability of a node having degree k 
% rho: threshold value
% 
% OUTPUTS
% alpha: power-law multiplier
% gamma: power-law exponent
% eps: mean fitting error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot of the probability distribution function
figure 
x = (0:max(deg));
plot(x,y,'-o');   
title(['Probability distribution with \rho = ', num2str(rho)]);

%% Least square fitting
[c,~,~] = lsqcurvefit_approx(x,y,'exp'); % suppose the distribution is a power-law 
                                         % with exponential fit  
xFit = linspace(min(x),max(x),50);
yFit = c(1)*exp(c(2)*xFit);

%% Plot the fitted probability distribution
hold on
plot(xFit,yFit,'-');
xlabel('k')
ylabel('p(k)')
legend('data','fitted');

%% Compute the parameters of the fitted probability distribution
alpha = c(1);
gamma = -c(2);

%% Compute the mean fitting error
error = y'-c(1)*exp(c(2)*x);
eps = sum(abs(error));       

%% Plot of the probability distribution function and the fitted one in loglog scale
figure
loglog(x,y,'-o');
xFit2 = linspace(1,max(x),50);
yFit2 = c(1)*exp(c(2)*xFit2);
hold on
loglog(xFit2,yFit2,'r-');
title (['Probability distribution with both axes in logarithmic scale and \rho = ', ...
        num2str(rho)]);
xlabel('k')
ylabel('p(k)')
legend('data','fitted');
end