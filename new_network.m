function [alpha,gamma,eps,new] = new_network(Corr,rho,names)
% function that creates the subnetwork and computes some basic quantities
% 
% INPUTS
% Corr: Pearson correlation matrix
% rho: threshold value
% names: list of the names of the stocks
% 
% OUTPUTS
% alpha: power-law multiplier
% gamma: power-law exponent
% eps: mean fitting error
% new: list of the positions of the stocks used to construct the subnetwork
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Construct the starting network 
dim = length(names);

A = adj_matrix(rho,Corr);     % adjacency matrix
B = abs(eye(dim)-1);
Adj = A.*B;     % elimination of selfloops
G = graph(Adj, names);
    
%% Compute some basic quantities of the subnetwork
[~,~,deg,prob,~,~,~,~,~,~,~,~,~] = graph_info(G,Adj,rho);
    
%% Compute the fitting model and the corresponding error 
[alpha,gamma,eps] = lsfit(deg,prob,rho);
    
%% Graph representation of the new subnetwork
new = find(deg);   % select the nodes with at least one neighbor
subG =  subgraph(G,new);
figure
plot(subG, 'NodeColor', 'k', 'EdgeColor', '[0 0.4470 0.7410]');
title(['Network with \rho = ', num2str(rho)]);
axis off
end