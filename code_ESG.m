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
load ESG.mat         % matrix with the ESG scores of each stock (column) for each time (row)

%% Sort the companies in decreasing order for ESG score 
sort_ESG = sort(mean(ESG),'descend');
sort_names_ESG = cell(1,length(sort_ESG));
for i = 1:length(sort_ESG)
    idx = find(mean(ESG) == sort_ESG(i));
    sort_names_ESG(i) = [names(idx)];
end

%% Correlation matrix analysis
[Indicator,max_ESG,min_ESG] = composite_indicator(Prices,ESG);
Corr_ESG = corrcoef(Indicator);    % Pearson correlation matrix
[max_negative_ESG,no_negative_ESG,max_value_ESG,high_corr_ESG] ...
                                  = corr_analysis(Corr_ESG,names); % other information

%% Network construction and parameters analysis
rho_ESG = [0, 0.3383, 0.5771, 0.8557, 0.9353, 0.9751];  % only interested in these values
                                                        % of the threshold vector 
dim = length(names);
N_ESG = zeros(length(rho_ESG),1);
L_ESG = zeros(length(rho_ESG),1);
deg_ESG = zeros(dim,length(rho_ESG));
d_ESG = zeros(length(rho_ESG),1);
k_ESG = zeros(length(rho_ESG),1);
h_ESG = zeros(length(rho_ESG),1);
D_ESG = zeros(length(rho_ESG),1);
avg_l_ESG = zeros(length(rho_ESG),1);
g_ESG = zeros(dim,length(rho_ESG));
b_ESG = zeros(dim,length(rho_ESG));
C_ESG = zeros(length(rho_ESG),1);
clustering_ESG = zeros(dim,length(rho_ESG));
alpha_ESG = zeros(length(rho_ESG),1);
gamma_ESG = zeros(length(rho_ESG),1);
eps_ESG = zeros(length(rho_ESG),1);
 
for i = 1:length(rho_ESG)
    A_ESG = adj_matrix(rho_ESG(i),Corr_ESG);      % adjacency matrix
    B_ESG = abs(eye(dim)-1);
    Adj_ESG = A_ESG.*B_ESG;    % elimination of selfloops
    G_ESG = graph(Adj_ESG, names);
    figure
    plot(G_ESG, 'NodeColor', 'k', 'EdgeColor', '[0 0.4470 0.7410]');
    title(['Network with \rho = ', num2str(rho_ESG(i))]);
    axis off
    
    % Compute some basic quantities of the network
    [N_ESG(i),L_ESG(i),deg_ESG(:,i),prob_ESG,d_ESG(i),k_ESG(i),h_ESG(i),D_ESG(i), ...
        avg_l_ESG(i),g_ESG(:,i),b_ESG(:,i),C_ESG(i), clustering_ESG(:,i)] ...
        = graph_info(G_ESG,Adj_ESG,rho_ESG(i));
    
    % Compute the fitting model and the corresponding error 
    [alpha_ESG(i),gamma_ESG(i),eps_ESG(i)] = lsfit(deg_ESG(:,i),prob_ESG,rho_ESG(i));
end

%% Plot and compute other informations about nodes centrality
% Count the number of stocks having the maximum degree for each example
[max_degree_ESG,n_max_degree_ESG] = centralities_info(deg_ESG,names);

% Count the number of stocks having the maximum closeness centrality for each example
[max_closeness_ESG,n_max_closeness_ESG] = centralities_info(g_ESG,names);

% Count the number of stocks having the maximum betweenness for each example
[max_betweenness_ESG,n_max_betweenness_ESG] = centralities_info(b_ESG,names);

% Count the number of stocks having the maximum clustering coefficient for each example
[max_clustering_ESG,n_max_clustering_ESG] = centralities_info(clustering_ESG,names);

%% Network construction in the scale-free range  
rho_range_ESG = linspace(0.9417,0.9664,50);    % threshold vector for scale-free network

for i = 1:length(rho_range_ESG)
    A_range_ESG = adj_matrix(rho_range_ESG(i),Corr_ESG);      % adjacency matrix
    B_range_ESG = abs(eye(dim)-1);
    Adj_range_ESG = A_range_ESG.*B_range_ESG;    % elimination of selfloops
    G_range_ESG = graph(Adj_range_ESG, names);
    
    % Compute some basic quantities of the network
    [~,~,deg_range_ESG,prob_range_ESG,~,~,~,~,~,~,~,~, ...
        ~] = graph_info(G_range_ESG,Adj_range_ESG,rho_range_ESG(i));
    
    % Compute the fitting model and the corresponding error 
    [~,~,~] = lsfit(deg_range_ESG,prob_range_ESG,rho_range_ESG(i));                                                   
end

%% Construct the new index  
Rho_ESG = [0.942848979591837 0.9542];    % same values as before 

alpha2_ESG = zeros(length(Rho_ESG),1);
gamma2_ESG = zeros(length(Rho_ESG),1);
eps2_ESG = zeros(length(Rho_ESG),1);
abs_err_new_50_ESG = zeros(size(Prices,1),length(Rho_ESG));   % su ogni colonna ho tt tempi
rel_err_new_50_ESG = zeros(size(Prices,1),length(Rho_ESG));

% Case rho = 0.9428
[alpha2_ESG(1),gamma2_ESG(1),eps2_ESG(1),new] = new_network(Corr_ESG,Rho_ESG(1),names);
[sort_names1_ESG,abs_err_new_50_ESG(:,1),rel_err_new_50_ESG(:,1)] =  ...
    new_index(new,Rho_ESG(1),Prices,Shares,Euro_price,dates,names);

% Case rho = 0.9542
[alpha2_ESG(2),gamma2_ESG(2),eps2_ESG(2),new] = new_network(Corr_ESG,Rho_ESG(2),names);
[sort_names2_ESG,abs_err_new_50_ESG(:,2),rel_err_new_50_ESG(:,2)] = ...
    new_index(new,Rho_ESG(2),Prices,Shares,Euro_price,dates,names);

%% Analysis of the stocks present in the new subnetworks
% Sort all the fifty companies in decreasing order for number of outstanding shares
sort_names_50_ESG = new_index((1:50),0,Prices,Shares,Euro_price,dates,names);
%%
% Plot the price evolution of the stocks remained when the index is
% approximated with few stocks for rho = 0.9428
price_evolution([31,27,36,48,32,35,8,13,6,21,20,41 ...
                 15,47,39,34,29,3,12,30,16,2,49],dates,Prices,names);
%%            
% Plot the price evolution of the stocks remained when the index is
% approximated with few stocks for rho = 0.9542
price_evolution([31,48,32,35,13,6,21,20,47,3,12,30,2,49],dates,Prices,names);
                                      
%% Error analysis
mean_abs_err_new_50_ESG = zeros(1,length(Rho_ESG));
var_abs_err_new_50_ESG = zeros(1,length(Rho_ESG));
mean_rel_err_new_50_ESG = zeros(1,length(Rho_ESG));
var_rel_err_new_50_ESG = zeros(1,length(Rho_ESG));

for i = 1:length(Rho_ESG)
    [mean_abs_err_new_50_ESG(i),var_abs_err_new_50_ESG(i)] = error_analysis(abs_err_new_50_ESG(:,i));
    [mean_rel_err_new_50_ESG(i),var_rel_err_new_50_ESG(i)] = error_analysis(rel_err_new_50_ESG(:,i));
end

%% Analysis of the stocks which differ between the case without and with the ESG score 
% Plot the price evolution of the stocks without the ESG score in case 1 which
% are not present in the same case with the ESG score
price_evolution([46,25,11,28,5,1,42],dates,Prices,names);
%%
% Plot the price evolution of the stocks with the ESG score in case 1 which
% are not present in the same case without the ESG score
price_evolution([36,13,41,39,29],dates,Prices,names);
%%
% Plot the price evolution of the stocks without the ESG score in case 2 which
% are not present in the same case with the ESG score
price_evolution([27,8,15,28,5,1],dates,Prices,names);
%%
% Plot the price evolution of the stocks with the ESG score in case 2 which
% are not present in the same case without the ESG score
price_evolution([13,6,30],dates,Prices,names);