function [deg,freq] = degree_distribution(G,rho)
% function that computes the degree of each node and plots the probability
% distribution function 
% 
% INPUTS
% G:  graph
% rho: threshold value
% 
% OUTPUTS
% deg: vector with the degree for each node
% freq: vector with the probability of a node having degree k
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deg = degree(G);     
count = zeros(max(deg)+1,1);    % number of nodes for each degree
for i = 1:length(deg)
    val = deg(i);
    count(val+1) = count(val+1)+1;  % in reality in position 1 there is degree 0
end
        
%% Degree distribution
s = sum(count);
freq = count/s;

%% Plot of the probability distribution function
figure 
x = (0:max(deg));
plot(x,freq,'-o');
title (['Probability distribution with \rho = ', num2str(rho)]);   
xlabel('k');
ylabel('p(k)');
end