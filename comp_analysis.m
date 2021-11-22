function [D,avg_l] = comp_analysis(G,comp)
% function that computes the diameter and the average path length for a
% disconnected network
% 
% INPUTS
% G: graph
% comp: subdivision of the network in components 
% 
% OUTPUTS
% D: diameter
% avg_l: average path length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diam = zeros(1,max(comp));       % vector with the diameter for each component
avg_leng = zeros(1,max(comp));   % vector with the average path length for each component
for i = 1:max(comp)
    new = find(comp==i);
    subG =  subgraph(G,new);
    N = numnodes(subG);   
    if N == 1
        diam(i) = 0;
        avg_leng(i) = 0;
    else
        L_max = N*(N-1)/2;
        l = distances(subG); % matrix with the length of the shortest path between 
                             % each pair of nodes
        diam(i) = max(l(:)); 
        T = triu(l);  
        avg_leng(i) = sum(T(:))/L_max;
    end
end
D = max(diam);
avg_l = mean(avg_leng(avg_leng~=0)); % consider only components with more than one element
end