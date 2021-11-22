function [N,L,deg,prob,d,k,h,D,avg_l,g,b,C,clustering] = graph_info(G,A,rho)
% function that computes some basic quantities of a graph G
% 
% INPUTS
% G: graph
% A: adjacency matrix
% rho: threshold value
% 
% OUTPUTS
% N: number of nodes
% L: number of links
% deg: vector with the degree for each node
% prob: vector with the probability of a node having degree k
% d: density 
% k: average degree 
% h: heterogeneity parameter
% D: diameter 
% avg_l: average path length
% g: vector with the closeness centrality for each node
% b: vector with the betweenness for each node
% C: average clustering coefficient
% clustering: vector with the clustering coefficient for each node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Basic quantities
N = numnodes(G);     
L = numedges(G);     
[deg,prob] = degree_distribution(G,rho);
L_max = N*(N-1)/2;
d = L/L_max;
k = d*(N-1);
h = mean(deg.^2)/k^2;

%% Component analysis
comp = conncomp(G); % check if there are different components
if max(comp)== 1   % if there is only a component
    l = distances(G);  % matrix with the length of the shortest path between
                       % each pair of nodes
    D = max(l(:));
    T = triu(l);    % upper triangular matrix
    avg_l = sum(T(:))/L_max;
else
    [D,avg_l] = comp_analysis(G,comp);      
end
    
%% Nodes quantities
g = centrality(G,'closeness');
b = centrality(G,'betweenness');
[C,clustering] = avgClusteringCoefficient(A);
end