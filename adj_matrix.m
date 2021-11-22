function [A] = adj_matrix(rho,C)
% function that computes the adjacency matrix 
% 
% INPUTS
% rho:  threshold used to select the links
% C: correlation matrix
% 
% OUTPUT
% A: adjacency matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_stocks = size(C,1);
A = zeros(n_stocks);
for i = 1:n_stocks
    for j = 1:n_stocks
        if abs(C(i,j)) > rho
            A(i,j) = 1;
        end
    end
end
end