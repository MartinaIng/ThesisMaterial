function [Indicator,max_ESG,min_ESG] = composite_indicator(Prices,ESG)
% function that computes the composite indicator 
%
% INPUTS
% Prices: matrix with the prices of each stock (column) for each time (row)    
% ESG: matrix with the ESG scores of each stock (column) for each time (row)
%
% OUTPUTS
% Indicator: matrix with the composite indicator of each stock (column) for each time (row)
% max_ESG: max ESG score
% min_ESG: min ESG score
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Normalization of the indicators
max_price = max(Prices);   % max price for each stock
max_ESG = max(ESG(:));     % max ESG score 
min_ESG = min(ESG(:));     % min ESG score

norm_Prices = zeros(size(Prices));
for i = 1:size(Prices,1)
    norm_Prices(i,:) = Prices(i,:)./max_price;  % standardized prices
end

norm_ESG = 1 - (max_ESG-ESG)/(max_ESG-min_ESG)*0.9;   % interpolation

%% Creation of the composite indicator
Indicator = sqrt(norm_Prices.*norm_ESG);   % geometric mean
end