function price_evolution(stocks,dates,Prices,names)
% function that plots the price evolution for different companies
%
% INPUTS
% stocks: vector of the positions of the stocks that are plotted
% dates: exchange dates
% Prices: matrix with the prices of each stock (column) for each time (row)    
% names: list of the names of the stocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = stocks
    figure
    plot(datetime(dates),Prices(:,i)/Prices(end,i))
    title(names(i));
    xlabel('time')
    ylabel('price')
end
end