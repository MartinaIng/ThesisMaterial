function [max_negative,no_negative,max_value,high_corr] = corr_analysis(Corr,names)
% function that plots the correlation matrix and computes some interesting 
% quantities related to correlation values
% 
% INPUTS
% Corr: Pearson correlation matrix
% names: list of the names of the stocks
%
% OUTPUTS
% max_negative: names of the stocks with the maximum number of negative 
%               correlaions followed by the value itself
% no_negative: names of the stocks with no negative correlation
% max_value: maximum value of the correlation matrix (apart from 1)
% high_corr: percentage of correlation coefficients with a value greater than 0.5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot of the correlation matrix
figure
imagesc(abs(Corr));   % display correlation matrix as an image
tickvalues = 1:length(Corr);
set(gca, 'YTick', tickvalues); 
set(gca, 'XTick', tickvalues);
x = zeros(size(tickvalues));
x(:) = length(Corr)+1;
text(tickvalues, x, names, 'HorizontalAlignment', 'right','Rotation',90); % set x-axis labels
set(gca, 'YTickLabel', names); % set y-axis labels
title('Correlation matrix');
colormap('jet');
colorbar; 

%% Count the number of negative correlation values for each stock
negative_corr = zeros(1,length(Corr));
for i = 1:length(Corr)
    col = Corr(:,i);
    negative_corr(i) = length(find(col < 0));
end

%% Names of the stocks with the maximum number of negative correlaions
M = max(negative_corr);    
idx = (negative_corr == M);
max_negative = [names(idx) M]; 

%% Names of the stocks with no negative correlation
P = (~negative_corr);   
no_negative = names(P);

%% Compute the maximum value of the correlation matrix (apart from 1)
corr = Corr-eye(length(Corr));
max_value = max(corr(:));

%% Percentage of the correlation coefficients over 0.5
high_corr = length(find(Corr > 0.5))*100/(size(Corr,1)*size(Corr,2));
end