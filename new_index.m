function [sort_names,abs_err_new_50,rel_err_new_50] = ...
          new_index(new,rho,Prices,Shares,Euro_price,dates,names)
% function that computes the new index and plots it against other prices
% along with the related absolute error, it computes also the absolute and
% the relative error between the new indices with 50 and less stocks  
%
% INPUTS
% new: list of the positions of the stocks used to construct the subnetwork
% rho: threshold value
% Prices: matrix with the prices of each stock (column) for each time (row)    
% Shares: matrix with the number of outstanding shares of each stock (column) 
%         for each time (row)
% Euro_price: vector with the prices of the Euro Stoxx 50 for each time 
% dates: exchange dates
% names: list of the names of the stocks
% 
% OUTPUTS
% sort_names: list of names of the stocks in decreasing order for number of 
%             outstanding shares 
% abs_err_new_50: absolute error between the new indeces with 50 and less 
%                 stocks prices 
% rel_err_new_50: relative error between the new indeces with 50 and less 
%                 stocks prices 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select prices and number of outstanding shares only for selected stocks 
Prices_subG = zeros(size(Prices,1),length(new));
Shares_subG = zeros(size(Shares,1),length(new));

i = 1; 
while i < length(new)+1
    Prices_subG(:,i) = Prices(:,new(i));
    Shares_subG(:,i) = Shares(:,new(i));
    names_subG(i) = names(new(i));
    i = i+1;
end

%% Sort the companies in decreasing order for number of outstanding shares 
sort_shares = sort(Shares_subG(end,:),'descend');
sort_names = cell(1,length(sort_shares));
for i = 1:length(sort_shares)
    idx = find(Shares_subG(end,:) == sort_shares(i));
    sort_names(i) = [names_subG(idx)];
end

%% New index computation
% with the stocks of the subnetwork
Index_new = sum(Prices_subG.*Shares_subG,2)/sum(Prices_subG(end,:).*Shares_subG(end,:),2);

% with all the 50 stocks
Index_new_50 = sum(Prices.*Shares,2)/sum(Prices(end,:).*Shares(end,:),2);

%% Plot the Euro Stoxx 50 and the new index prices evolutions (with less stocks)
figure
plot(datetime(dates),Euro_price/Euro_price(end));
xlabel('time')
ylabel('price')
title(['Comparison between the two indices with \rho = ', num2str(rho)]);
hold on
plot(datetime(dates),Index_new);
legend('Euro Stoxx 50','New index','Location','northwest');

% Plot of the absolute error
abs_err_new = abs(Euro_price/Euro_price(end)-Index_new); 
figure
plot(datetime(dates),abs_err_new);
xlabel('time')
ylabel('error')
title(['Error with \rho = ', num2str(rho)]);

%% Plot both the new indeces prices evolutions (with 50 and less stocks)
figure
plot(datetime(dates),Index_new_50);
xlabel('time')
ylabel('price')
title(['Comparison between the two indices with \rho = ', num2str(rho)]);
hold on
plot(datetime(dates),Index_new);
legend('Index with 50 stocks','Index with less stocks','Location','northwest');

% Plot of the absolute error
abs_err_new_50 = abs(Index_new_50-Index_new);
rel_err_new_50 = abs_err_new_50./Index_new_50;
figure
plot(datetime(dates),abs_err_new_50);
xlabel('time')
ylabel('error')
title(['Error with \rho = ', num2str(rho)]);
end