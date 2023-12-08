%% Compute the SADF and GSADF statistics for the price-dividend ratio.
% This file contains a sample run for the procedure developed in:

%Tomas E. Caravello, Zacharias Psaradakis, Martin Sola,
% "Rational bubbles: Too many to be true?", Journal of Economic Dynamics and Control, Volume 151, 2023
%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/tomyc/Documents/GitHub/bubbles_bootstrap/code'; % replace here with your path.


addpath([path '/_aux'])
addpath([path '/_data'])

cd(path);
%% Load stuff
% import
SPDV= xlsread('shiller_quarterly.xlsx','Sheet1'); %Shiller data on stock prices and dividends, quarterly.
%%
year = SPDV(226:596,1);
y = SPDV(226:596,2)./SPDV(226:596,3); %data series to use. In this case, we use the price dividend ratio for a subset of observations.

% column 1: year
% column 2: real stock prices
% column 3: real dividends (annualized)

T=length(y); 
r0=0.01+1.8/sqrt(T);  
swindow0 = floor(r0*T);
dim=T-swindow0+1;
nrep_b = 10000; % nreb_b is used for the BSADF sequence. Should be moderately large to achieve good precision at the tails.
nrep2 = 1000;  % number of iters for the first stage of the stwo step procedure. Doesn't need to be that large. 
M = 500; % number of repetitions for the second stage in the bootstrapping. 
swindow0 = floor(r0*T);
dim=T-swindow0+1;
paralell = 1;
verbose_b = 0;

quantile_sadf_asymp = [1.20;1.49; 2.07]; %quantiles obtained from PSY (2011) 
quantile_gsadf_asymp = [1.97; 2.19; 2.69];
%%

% SEE DOCUMENATION OF THE FUNCTION. Many of the outputs we used for
% simulations, but may be not relevant for applied use.  

[sadf,badfs,gsadf,bsadfs,...
            cvs_sadf_calibrated,cvs_gsadf_calibrated,...
            criticalvs_badf,criticalvs_bsadf,...
           reject_sadf_asym,reject_gsadf_asym,reject_sadf_not_calibrated,...
           reject_gsadf_not_calibrated,reject_sadf_calibrated,reject_gsadf_calibrated,...
           reject_badf_no_calib_no_precheck, reject_bsadf_no_calib_no_precheck,...  
           reject_badf_calib_no_precheck, reject_bsadf_calib_no_precheck,...
           reject_badf_no_calib_precheck, reject_bsadf_no_calib_precheck,...  
           reject_badf_calib_precheck, reject_bsadf_calib_precheck] =...
    run_everything(y,r0,nrep_b,nrep2,M,paralell,verbose_b, 12, 0.05,quantile_sadf_asymp,quantile_gsadf_asymp,1);

%% Print results and plot
clc

fprintf(' Statistic & \\multicolumn{3}{c}{Asymptotic} & & \\multicolumn{3}{c}{Calibrated} \\\\ \n')
fprintf(' &  & 0.1 & 0.05 & 0.01 & & 0.1 & 0.05 & 0.01 \\\\ \n')
fprintf('SADF & %2.4f  & %2.4f & %2.4f & %2.4f && %2.4f & %2.4f & %2.4f \\\\ \n', [sadf; quantile_sadf_asymp ; cvs_sadf_calibrated])
fprintf('GSADF & %2.4f  & %2.4f & %2.4f & %2.4f && %2.4f & %2.4f & %2.4f \n', [gsadf;quantile_gsadf_asymp ; cvs_gsadf_calibrated])


%% Plot
close all
figure (1)
plot(year(swindow0:end),bsadfs','LineWidth',4)
hold on
plot(year(swindow0:end),criticalvs_bsadf,'LineWidth',4);
hold off
%axis([1929 2020 -1.5 4.5])
%xticks(1929:10:2020)
title('The backward SADF sequence and 95% critical values - Price-dividend ratio','FontSize',10);
xlabel('Year');
ylabel('BSADF');
legend('Statistic','Critical value - 95%','location','northwest','orientation','vertical')
