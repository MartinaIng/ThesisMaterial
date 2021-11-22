clear all;
close all;
clc;

%% Upload the dataset (dates are in decreasing order)
load dates.mat          % exchange dates
load names.mat          % names of the stocks

load Prices.mat      % matrix with the prices of each stock (column) for each time (row)  
load Euro_price.mat  % vector with the prices of the Euro Stoxx 50 for each time 
load Shares.mat      % matrix with the number of outstanding shares of each stock (column) 
                     % for each time (row)

%% Correlation matrix analysis
Corr = corrcoef(Prices);    % Pearson correlation matrix
[max_negative,no_negative,max_value,high_corr] = corr_analysis(Corr,names); % other information

%% Network construction and parameters analysis
rho = linspace(0,0.9751,50);    % threshold vector   

dim = length(names);
N = zeros(length(rho),1);
L = zeros(length(rho),1);
deg = zeros(dim,length(rho));
d = zeros(length(rho),1);
k = zeros(length(rho),1);
h = zeros(length(rho),1);
D = zeros(length(rho),1);
avg_l = zeros(length(rho),1);
g = zeros(dim,length(rho));
b = zeros(dim,length(rho));
C = zeros(length(rho),1);
clustering = zeros(dim,length(rho));
alpha = zeros(length(rho),1);
gamma = zeros(length(rho),1);
eps = zeros(length(rho),1);

for i = 1:length(rho)
    A = adj_matrix(rho(i),Corr);      % adjacency matrix
    B = abs(eye(dim)-1);
    Adj = A.*B;    % elimination of selfloops
    G = graph(Adj, names);
    figure
    plot(G, 'NodeColor', 'k', 'EdgeColor', '[0 0.4470 0.7410]');
    title(['Network with \rho = ', num2str(rho(i))]);
    axis off
    
    % Compute some basic quantities of the network
    [N(i),L(i),deg(:,i),prob,d(i),k(i),h(i),D(i),avg_l(i),g(:,i),b(:,i),C(i), ...
        clustering(:,i)] = graph_info(G,Adj,rho(i));
    
    % Compute the fitting model and the corresponding error 
    [alpha(i),gamma(i),eps(i)] = lsfit(deg(:,i),prob,rho(i));
end

%% Plot and compute other informations about nodes centrality
% Count the number of stocks having the maximum degree for each example
[max_degree,n_max_degree] = centralities_info(deg,names);

% Count the number of stocks having the maximum closeness centrality for each example
[max_closeness,n_max_closeness] = centralities_info(g,names);

% Count the number of stocks having the maximum betweenness for each example
[max_betweenness,n_max_betweenness] = centralities_info(b,names);

% Count the number of stocks having the maximum clustering coefficient for each example
[max_clustering,n_max_clustering] = centralities_info(clustering,names);

%% Network construction in the scale-free range
rho_range = linspace(0.9308,0.9636,50);    % threshold vector for scale-free network

h_range = zeros(length(rho_range),1);
eps_range = zeros(length(rho_range),1);
LCC = zeros(length(rho_range),1);
Eff = zeros(length(rho_range),1);

for i = 1:length(rho_range)
    A_range = adj_matrix(rho_range(i),Corr);      % adjacency matrix
    B_range = abs(eye(dim)-1);
    Adj_range = A_range.*B_range;    % elimination of selfloops
    G_range = graph(Adj_range, names);
    
    % Compute some basic quantities of the network
    [~,~,deg_range,prob_range,~,~,h_range(i),~,~,~,~,~, ...
        ~] = graph_info(G_range,Adj_range,rho_range(i));
    
    % Compute the fitting model and the corresponding error 
    [~,~,eps_range(i)] = lsfit(deg_range,prob_range,rho_range(i));
                                                         
    % Compute the Largest Connected Cluster and the Efficiency functions
    [LCC(i),Eff(i)] = LCC_Eff(G_range);
end

%% Plot the heterogeneity parameter evolution in the scale-free range
figure
plot(rho_range,h_range,'-o');
title('Heterogeneity evolution');
xlabel('\rho');
ylabel('h');

%% Plot the mean fitting error in the scale-free range
figure
plot(rho_range,eps_range,'-o');
title('Mean fitting error evolution');
xlabel('\rho');
ylabel([char(949) '_{fitting}']);

%% Plot LCC and Efficiency evolution in the scale-free range
figure
plot(rho_range,LCC,'-o');
title('Largest Connected Cluster evolution');
xlabel('\rho');
ylabel('LCC');

figure
plot(rho_range,Eff,'-o');
title('Efficiency evolution');
xlabel('\rho');
ylabel('Eff');

%% Construct the new index
idx = find(eps_range == min(eps_range));    % index of the rho with less error
Rho = [0.9355 rho_range(idx) 0.9542];       % considered rho for subgraph construction

alpha2 = zeros(length(Rho),1);
gamma2 = zeros(length(Rho),1);
eps2 = zeros(length(Rho),1);
abs_err_new_50 = zeros(size(Prices,1),length(Rho));   
rel_err_new_50 = zeros(size(Prices,1),length(Rho)); 

% Case rho = 0.9355
[alpha2(1),gamma2(1),eps2(1),new] = new_network(Corr,Rho(1),names);
[sort_names1,abs_err_new_50(:,1),rel_err_new_50(:,1)] = ...
    new_index(new,Rho(1),Prices,Shares,Euro_price,dates,names);

% Case rho = 0.9428
[alpha2(2),gamma2(2),eps2(2),new] = new_network(Corr,Rho(2),names);
[sort_names2,abs_err_new_50(:,2),rel_err_new_50(:,2)] = ...
    new_index(new,Rho(2),Prices,Shares,Euro_price,dates,names);

% Case rho = 0.9542
[alpha2(3),gamma2(3),eps2(3),new] = new_network(Corr,Rho(3),names);
[sort_names3,abs_err_new_50(:,3),rel_err_new_50(:,3)] = ...
    new_index(new,Rho(3),Prices,Shares,Euro_price,dates,names);

%% Analysis of the deleted stocks
% Sort all the fifty companies in decreasing order for number of outstanding shares
sort_names_50 = new_index((1:50),0,Prices,Shares,Euro_price,dates,names);
%%
% Plot the price evolution of the stocks eliminated from case Euro Stoxx 50 to case 1
price_evolution([19,24,44,14,10,50,43,13,9,4,41,40,33,23, ...
                 37,39,45,22,18,26,38,17],dates,Prices,names);
%%
% Plot the price evolution of the stocks eliminated from case 1 to case 2
price_evolution([36,29,7],dates,Prices,names);
%%
% Plot the price evolution of the stocks eliminated from case 2 to case 3
price_evolution([46,25,6,11,34,42,30,16],dates,Prices,names);

%% Error analysis 
mean_abs_err_new_50 = zeros(1,length(Rho));
var_abs_err_new_50 = zeros(1,length(Rho));
mean_rel_err_new_50 = zeros(1,length(Rho));
var_rel_err_new_50 = zeros(1,length(Rho));

for i = 1:length(Rho)    
    [mean_abs_err_new_50(i),var_abs_err_new_50(i)] = error_analysis(abs_err_new_50(:,i));
    [mean_rel_err_new_50(i),var_rel_err_new_50(i)] = error_analysis(rel_err_new_50(:,i));
end