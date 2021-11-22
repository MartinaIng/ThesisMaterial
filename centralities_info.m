function[max_value,n_max] = centralities_info(matrix,names)
% function that plots, for each example, the centrality measure of each node 
% and that computes some interesting quantities related to the maximum values reached
% 
% INPUTS
% matrix: matrix with the quantity for each node (rows) for each example (columns)
% names: list of the names of the stocks
%
% OUTPUTS
% max_value: vector with the maximum value of the analyzed quantity for each example
% n_max: number of nodes that have the highest value 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_value = max(matrix); 
n_max = zeros(1,length(max_value));
figure
x = 1:length(names);
for i = 1:length(max_value)
    col = matrix(:,i);
    idx = find(col == max_value(i));
    n_max(i) = length(idx);   
    grid on
    hold on
    plot(x,col,'LineWidth',1.5);
end
tickvalues = x;
set(gca, 'XTick', tickvalues); 
x = zeros(size(tickvalues));
text(tickvalues, x, names, 'HorizontalAlignment', 'right','Rotation',90);
end