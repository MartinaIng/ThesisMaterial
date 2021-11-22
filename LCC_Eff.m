function [LCC, Eff] = LCC_Eff(G)
% function that computes the Largest Connected Cluster and the Efficiency 
% of a graph and plot them
% 
% INPUT
% G: graph
%
% OUTPUTS
% LCC: maximum number of connected nodes
% Eff: efficiency coefficient of the graph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute the Largest Connected Cluster 
comp = conncomp(G);     % divides the nodes in different components
count = zeros(max(comp),1);    % number of nodes for each component
for i = 1:length(comp)
    val = comp(i);
    count(val) = count(val)+1;  
end
LCC = max(count);

%% Compute the Efficiency 
N = numnodes(G); 
l = diag(Inf*ones(1,N)) + distances(G);   % matrix with the length of the 
                                          % shortest path between each pair
                                          % of nodes, Inf for a node with itself 
l_inv = 1./l;   
T = triu(l_inv);         
Eff = sum(T(:))/(N*(N-1));
end