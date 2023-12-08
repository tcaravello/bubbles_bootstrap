function  [sadf,badfs,gsadf,bsadfs,...
            cvs_sadf_calibrated,cvs_gsadf_calibrated,...
            criticalvs_badf,criticalvs_bsadf,...
    reject_sadf_asym,reject_gsadf_asym,reject_sadf_not_calibrated,...
           reject_gsadf_not_calibrated,reject_sadf_calibrated,reject_gsadf_calibrated,...
           reject_badf_no_calib_no_precheck, reject_bsadf_no_calib_no_precheck,...  
           reject_badf_calib_no_precheck, reject_bsadf_calib_no_precheck,...
           reject_badf_no_calib_precheck, reject_bsadf_no_calib_precheck,...  
           reject_badf_calib_precheck, reject_bsadf_calib_precheck,...
            lambdas_sadf, lambdas_gsadf,...
            calibrated_significance_level_badf,calibrated_significance_level_bsadf] =...
    run_everything(y,r0,nrep_b,nrep2,M,paralell,verbose_b, k, significance_bsadf,quantile_sadf_asymp,quantile_gsadf_asymp,robust)

% INPUTS
% - data: T x 1 data vector
% - r0: minimum fraction fo sample to use
% - nrep_b: number of bootstrap repetitions, first pass. This one should be large, since it requires precision 
% nreb_b is used for the BSADF sequence
% - nrep2: % number of iters for the first stage of the stwo step procedure. Doesn't need to be that large.
% - M:  number of repetitions for the second stage in the bootstrapping.
% - paralell: 1 if you wanna paralelize the function, 0 otherwise. For
% Monte Carlo, put at 0 and do the paralelization in the outer loop.
% - verbose_b: if you wanna print intermediate iterations for the
% bootstrap, set at 1. Otherwise, at 0.
% - k: minimum number of observations above the threshold to be considered
% a bubble.
% - significance_bsadf: significance level for the bsadf. You probably
% should use 0.05
% - quantile_sadf_asymp: 3x1 vector of asymptotic CVs
% - quantile_gsadf_asymp. 3x1 vector of asymptotic CVs
% - robust = 1 if you want use heteorskedasticity robust bootstrap. 0
% otherwise

% OUTPUTS
% Here pre check means that, before date-stamping, we check if we rejected
% the null for the whole sample. That is, we check if GSADF> critical
% value. 

% sadf: 1 x 1, sadf in the data,
% badfs: dim x 1, badf sequence
% gsadf: 1 x 1, gsadf in the data,
% bsadfs: 1 x dim, bsadf sequence
% cvs_sadf_calibrated: 3 x 1, sadf double bootstrap calibrated critical values for sizes 10, 5 and 1 per cent
% cvs_gsadf_calibrated: 3 x 1, gsadf double bootstrap calibrated critical values for sizes 10, 5 and 1 per cent
% criticalvs_badf: dim x 1, critical values for badf, calibrated to account
% for multiple testing.
% criticalvs_bsadf: dim x 1, critical values for badf, calibrated to account
% for multiple testing.
% reject_sadf_asym: 3 x 1, logical. 1 indicates rejection with aympt CV
% reject_gsadf_asym:  3 x 1, logical. 1 indicates rejection with aympt CV
% reject_sadf_not_calibrated: 3 x 1, logical. 1 indicates rejection with bootstrapped CV, not calibrated
% reject_gsadf_not_calibrated: 3 x 1, logical. 1 indicates rejection with bootstrapped CV, not calibrated
% reject_sadf_calibrated: 3 x 1, logical. 1 indicates rejection with bootstrapped CV, calibrated with double bootstrap
% reject_gsadf_calibrated: 3 x 1, logical. 1 indicates rejection with bootstrapped CV, calibrated with double bootstrap
% reject_badf_no_calib_no_precheck: 1 x 1, number of explosive episodes detected. CVs: not calibrated for multiple testing, no precheck
% reject_bsadf_no_calib_no_precheck : 1 x 1, number of explosive episodes detected. CVs: not calibrated for multiple testing, no precheck
% reject_badf_calib_no_precheck: 1 x 1, number of explosive episodes detected. CVs: calibrated for multiple testing (i.e, using criticalvs_badf), no precheck
% reject_bsadf_calib_no_precheck:  1 x 1, number of explosive episodes detected. CVs: calibrated for multiple testing (i.e, using criticalvs_bsadf), no precheck
% reject_badf_no_calib_precheck:  1 x 1, number of explosive episodes detected. CVs: not calibrated for multiple testing, precheck with asymptotic critical values.
% reject_bsadf_no_calib_precheck:  1 x 1, number of explosive episodes detected. CVs: not calibrated for multiple testing, precheck with asymptotic critical values.
% reject_badf_calib_precheck: 1 x 1, number of explosive episodes detected.
% CVs: calibrated for multiple testing (i.e, using criticalvs_bsadf), precheck with double-boostrap critical values (i.e, with cvs_sadf_calibrated) 
% reject_bsadf_calib_precheck: 1 x 1, number of explosive episodes detected.
% CVs: calibrated for multiple testing (i.e, using criticalvs_bsadf), precheck with double-boostrap critical values (i.e, with cvs_sadf_calibrated) 
% lambdas_sadf: 3 x 1, calibrated significance level lambda in order to get a true rejection rate of 10,5 and 1 percent respectively. 
% lambdas_gsadf: 3 x 1, calibrated significance level lambda in order to get a true rejection rate of 10,5 and 1 percent respectively. 
% calibrated_significance_level_badf: 1 x 1, calibrated significance level lambda in order to get a true rejection rate of 5 per cent, accounting for multiple testing. 
% calibrated_significance_level_bsadf: 1 x 1, calibrated significance level lambda in order to get a true rejection rate of 5 per cent, accounting for multiple testing. 

% Written by Tomas Caravello, 12/20/2022

T = length(y);
swindow0=floor(r0*T);
dim=T-swindow0+1;

% First, get stuff from the actual data.

[sadf,badfs,gsadf,bsadfs] = get_SADF_GSADF_and_sequences(y,T,dim,swindow0);

% Two procedures to be implemented
% 1) correction for the trajectories of BADF/BSADF
% 2) correction for the smaple statistic SADF/GSADF

% for 1) we need to bootstrap from the true sample. Run first bootstrap iteration and get valuable inputs
[empirsadf,empirgsadf,traject_bootstrap,traject_badf,traject_bsadf] = ...
    bootstrap_from_sample_adjusted_cvs(y,nrep_b,T,swindow0,dim,verbose_b,paralell,robust); 

%% with these objects, I can run both 1) and 2)  

% obtain critical cvs for trajectories from bootstrapping
[criticalvs_badf,criticalvs_bsadf,calibrated_significance_level_badf,calibrated_significance_level_bsadf] =...
    obtain_cvs_for_bootstrap(traject_badf,traject_bsadf,k,significance_bsadf);

% obtain critical cvs for trajectories from bootstrapping with pre check

qe = [0.90,0.95,0.99];
% obtain cvs from bootstrap without calibration.
quantile_sadf_no_calibration = (quantile(empirsadf',qe))';
quantile_gsadf_no_calibration = (quantile(empirgsadf',qe))';

% obtain trajectory of badfs and bsadfs without calibration
quantile_badfs_no_calib = (quantile(traject_badf',qe))';
quantile_bsadfs_no_calib = (quantile(traject_bsadf',qe))';

%% run 2)
% run the bootstrapping par
data_to_use = traject_bootstrap(1:end,1:nrep2);
[colector_sadf_distribution,colector_gsadf_distribution,colector_sadf_true_values,colector_gsadf_true_values] ...
    = run_double_bootstrap(data_to_use, M,T,swindow0,dim, paralell,verbose_b,robust);

% get calibrated critical values.
[cvs_sadf_calibrated,cvs_gsadf_calibrated,  lambdas_sadf, lambdas_gsadf] = double_bootstrap_analyse_results(colector_sadf_distribution,colector_gsadf_distribution,...
    colector_sadf_true_values,colector_gsadf_true_values,...
    empirgsadf,empirsadf);

% use those to compute the sequences of badf and bsadf with precheck.

%[quantile_badfs_calib_precheck,quantile_bsadfs_calib_precheck,significance_level_badf,significance_level_bsadf] ...
 %   = obtain_cvs_for_bootstrap_precheck(traject_badf,traject_bsadf,k,cvs_sadf_calibrated(2),cvs_gsadf_calibrated(2),significance_bsadf)



%using all computed objects, check whether we reject or not.
%% Check rejections

% whole sample statistics
reject_sadf_asym = sadf*ones(3,1)>quantile_sadf_asymp';
reject_gsadf_asym = gsadf*ones(3,1)>quantile_gsadf_asymp';
reject_sadf_not_calibrated = sadf*ones(3,1)>quantile_sadf_no_calibration;
reject_gsadf_not_calibrated = gsadf*ones(3,1)>quantile_gsadf_no_calibration;
reject_sadf_calibrated = gsadf*ones(3,1)>cvs_sadf_calibrated;
reject_gsadf_calibrated = gsadf*ones(3,1)>cvs_gsadf_calibrated;

% remember to output the lambdas

% badf/bsadf sequences. All these are the number of detected bubbles
% note: the pattern of transposes for badfs and bsadf comes just from how
% stuff are defined in the code.

% no precheck: don't look at rejection with the SADF/GSADF.
reject_badf_no_calib_no_precheck = compute_number_of_bubbles(badfs,quantile_badfs_no_calib(:,2),k);
reject_bsadf_no_calib_no_precheck = compute_number_of_bubbles(bsadfs',quantile_bsadfs_no_calib(:,2),k);
reject_badf_calib_no_precheck = compute_number_of_bubbles(badfs,criticalvs_badf,k);   
reject_bsadf_calib_no_precheck = compute_number_of_bubbles(bsadfs',criticalvs_bsadf,k);

% precheck: look at rejection for SADF/GSADF, with asymptotic critical values. 
% not calibrated
if reject_sadf_asym(2) == 0
    reject_badf_no_calib_precheck = 0;
else
    reject_badf_no_calib_precheck = reject_badf_no_calib_no_precheck;
end

if reject_gsadf_asym(2) == 0
    reject_bsadf_no_calib_precheck = 0;
else
    reject_bsadf_no_calib_precheck = reject_bsadf_no_calib_no_precheck;
end

% calibrated
if reject_sadf_calibrated(2) == 0
    reject_badf_calib_precheck = 0;
else
    reject_badf_calib_precheck = reject_badf_calib_no_precheck; 
end

if reject_gsadf_calibrated(2) == 0
    reject_bsadf_calib_precheck = 0;
else
    reject_bsadf_calib_precheck = reject_bsadf_calib_no_precheck; 
end

end

% All functions called by the main function are below for extra speed.

%% Functions
% this file computes the sequence of BSDAF in sample for an AR(2) model. 
function [sadf,badfs,gsadf,bsadfs] = get_SADF_GSADF_and_sequences(y,T,dim,swindow0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the SADF and GSADF statistic, and the
% corresponsing BADF and BSADF sequences, by running an AR(2) as suggested
% by Wang and Yu (2022).

% Inputs:
% - y: Tx1 data vector
% - r0: minium fraction of the sample required to run the regression.
% Recommended: r0 = 0.01+1.8/sqrt(T)
% Outputs
% - sadf: value of the SADF statistic
% - gsadf: value of the GSADF statistic
% - badf: sequence of backward ADFs, from the computatio of the SADF
% - bsadf: sequence of backward SADFs, from the computatio of the GSADF

% written by Tomas Caravello. 
% first version: 12/19/2022
% this version: 12/19/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Actual Data
% SADF test and backward ADFs statistics.

badfs = ADF_FL4(y,swindow0,T);
sadf= max(badfs);

% GSADF and BSADFs
r2=(swindow0:1:T)';
bsadfs=zeros(1,dim);
for v=1:1:size(r2,1);
    rwadft = ADF_FL4_backwards(y(1:r2(v)),swindow0,r2(v));
    bsadfs(1,v)=max(rwadft);
 end;
gsadf=max(bsadfs(1,:),[],2);
end



function estm_seq = ADF_FL4(y,swindow0,T)

    % let y be the full data vector, swindow0 the smalles sample size to
    % use. 
    %step 1: demean everything properly.
    %mean_y = cumsum(y(1:T))./((1:T)');
    % truncate the means
    %mean_y(1:swindow0-1) = mean_y(swindow0);
    mean_y0 = cumsum(y(2:T))./((1:T-1)');
    mean_y1 = cumsum(y(1:T-1))./((1:T-1)');
    mean_y0(1:swindow0-2) = mean_y0(swindow0-1);
    mean_y1(1:swindow0-2) = mean_y1(swindow0-1);
    
    %y0 = y(2:T)-mean_y0;
    %y1 = y(1:T-1)-mean_y1;
    
    %xx_mat = cumsum(y1.^2);
    %yy_mat = cumsum(y0.^2);
    %xy_mat = cumsum(y0.*y1);

    yy_mat = cumsum(y(2:T).^2)-(mean_y0.^2).*(1:T-1)' ;
    xx_mat = cumsum(y(1:T-1).^2)-(mean_y1.^2).*(1:T-1)';
    xy_mat = cumsum(y(1:T-1).*y(2:T))-(mean_y0.*mean_y1).*(1:T-1)';


  beta = xy_mat(swindow0-1:T-1)./xx_mat(swindow0-1:T-1);                       % model A-@
 
  sig = sqrt((yy_mat(swindow0-1:T-1)./xx_mat((swindow0-1:T-1))).*(1-((xy_mat(swindow0-1:T-1).^2)./(xx_mat(swindow0-1:T-1).*yy_mat(swindow0-1:T-1))))./((swindow0-3:T-3)'));
  
  %if mflag==1; sig = sqrt(diag(eps'*eps/(t2-adflag-2)*(x2'*x2)^(-1))); end;
  %if mflag==2; sig = sqrt(diag(eps'*eps/(t2-adflag-3)*(x2'*x2)^(-1))); end;
  %if mflag==3; sig = sqrt(diag(eps'*eps/(t2-adflag-1)*(x2'*x2)^(-1))); end;

  %tvalue=beta./sig;
 
  estm_seq=(beta-1)./sig;

end


function estm_seq = ADF_FL4_backwards(y,swindow0,T)

    % let y be the full data vector, swindow0 the smalles sample size to
    % use. 
    %step 1: demean everything properly.
    %mean_y = cumsum(y(1:T))./((1:T)');
    % truncate the means
    %mean_y(1:swindow0-1) = mean_y(swindow0);
    mean_y0 = cumsum(flip(y(2:T)))./((1:T-1)');
    mean_y1 = cumsum(flip(y(1:T-1)))./((1:T-1)');
    mean_y0(1:swindow0-2) = mean_y0(swindow0-1);
    mean_y1(1:swindow0-2) = mean_y1(swindow0-1);
    
    %y0 = y(2:T)-mean_y0;
    %y1 = y(1:T-1)-mean_y1;
    
    %xx_mat = cumsum(y1.^2);
    %yy_mat = cumsum(y0.^2);
    %xy_mat = cumsum(y0.*y1);

    yy_mat = cumsum(flip(y(2:T)).^2)-(mean_y0.^2).*(1:T-1)' ;
    xx_mat = cumsum(flip(y(1:T-1)).^2)-(mean_y1.^2).*(1:T-1)';
    xy_mat = cumsum(flip(y(1:T-1)).*flip(y(2:T)))-(mean_y0.*mean_y1).*(1:T-1)';


  beta = xy_mat(swindow0-1:T-1)./xx_mat(swindow0-1:T-1);                       % model A-@
 
  sig = sqrt((yy_mat(swindow0-1:T-1)./xx_mat((swindow0-1:T-1))).*(1-((xy_mat(swindow0-1:T-1).^2)./(xx_mat(swindow0-1:T-1).*yy_mat(swindow0-1:T-1))))./((swindow0-3:T-3)'));
  
  %if mflag==1; sig = sqrt(diag(eps'*eps/(t2-adflag-2)*(x2'*x2)^(-1))); end;
  %if mflag==2; sig = sqrt(diag(eps'*eps/(t2-adflag-3)*(x2'*x2)^(-1))); end;
  %if mflag==3; sig = sqrt(diag(eps'*eps/(t2-adflag-1)*(x2'*x2)^(-1))); end;

  %tvalue=beta./sig;
 
  estm_seq=(beta-1)./sig;
end




function [empirsadf,empirgsadf,traject_bootstrap,traject_badf,traject_bsadf] = ...
    bootstrap_from_sample_adjusted_cvs(data,nrep_b,T,swindow0,dim,verbose_b,paralell,robust) 

% This function computes the bootstrapped distribution for the statistic
% for a given sample data or size T.
% INPUTS
% - data: series (in levels) to bootstrap. Size T
% - nrep_b: number of bootstrap iterations
% - T: sample size of data
% - swindow0: minium number of observations to run the procedure.
% - dim: dim=T-swindow0+1
% - verbose_b = 1 if you want to print time for each iteration

% OUTPUT
% - trayec_bootsrap: a matrix that contains the bootstrapped trayectories
% - empirsadf: the bootstrapped distribution for the SADF
% - empirgsadf: the bootstrapped distribution for the GSADF
% - true_sadf: the value of the SADF in the sample
% - true_gsadf: the value of the GSADF in the sample
%  - collect

% Create space to store some things
traject_bootstrap = zeros(T,nrep_b); % store bootstraped trajectories
traject_badf = zeros(dim,nrep_b); % store bootstraped BADF trajectories
traject_bsadf = zeros(dim,nrep_b); % store bootstraped BSADF trajectories
empirsadf = zeros(nrep_b,1); % store the SADF statistic
empirgsadf = zeros(nrep_b,1); % store the GSADF statistic 

data_init = data(1,1);
ret = data(2:end)-data(1:end-1);
alpha = mean(ret);
resid = ret-alpha;
sigma_hat = std(ret);
% INIT MAIN LOOP
if paralell == 1
   parfor nn_b = 1:nrep_b
        if verbose_b ==1
        tic
    fprintf('Stating bootstrap iter %f \n',nn_b)
        end
       % Step 1: create bootstraped trayectories as explained in the paper
          % Step 1: create bootstraped trayectories as explained in the paper
          if robust == 0
              deltay = alpha+sigma_hat*randn(size(resid));
          else
              deltay = alpha+resid.*randn(size(resid));
          end
        y = [data_init; data_init+cumsum(deltay)] ;
        traject_bootstrap(:,nn_b) = y;
    
       % Compute the statistics
        %THE SUP ADF TEST %%%%%%
        badfs= ADF_FL4(y,swindow0,T);
        sadf=max(badfs);
        traject_badf(:,nn_b) = badfs;
        empirsadf(nn_b,1)=sadf; 
        % badfs: vector with the ADFs. This trajectory will be compared with
        % the critical values.
    %THE GENERALIZED SUP ADF TEST %%%%%%
         
        r2=(swindow0:1:T)';
        bsadfs=zeros(1,dim);
        for v=1:1:size(r2,1);
            rwadft = ADF_FL4_backwards(y(1:r2(v)),swindow0,r2(v));
            bsadfs(1,v)=max(rwadft);
        end;
        gsadf=max(bsadfs(1,:),[],2);
        traject_bsadf(:,nn_b) = bsadfs'; 
        % bsadfs: vector with the SADFs. This trajectory will be compared with
        % the critical values.
        empirgsadf(nn_b,1) = gsadf;
    if verbose_b == 1;
        toc
        clc
    end
    
    end

else
    for nn_b = 1:nrep_b
        if verbose_b ==1
        tic
    fprintf('Stating bootstrap iter %f \n',nn_b)
        end
       % Step 1: create bootstraped trayectories as explained in the paper
          if robust ==0
              deltay = alpha+sigma_hat*randn(size(resid));
          else
              deltay = alpha+resid.*randn(size(resid));
          end
       y = [data_init; data_init+cumsum(deltay)] ; 
       traject_bootstrap(:,nn_b) = y;
    
       % Compute the statistics
        %THE SUP ADF TEST %%%%%%
        badfs= ADF_FL4(y,swindow0,T);
        sadf=max(badfs);
        traject_badf(:,nn_b) = badfs;
        empirsadf(nn_b,1)=sadf; 
        % badfs: vector with the ADFs. This trajectory will be compared with
        % the critical values.
    %THE GENERALIZED SUP ADF TEST %%%%%%
         
        r2=(swindow0:1:T)';
        bsadfs=zeros(1,dim);
        for v=1:1:size(r2,1);
            rwadft = ADF_FL4_backwards(y(1:r2(v)),swindow0,r2(v));
            bsadfs(1,v)=max(rwadft);
        end;
        gsadf=max(bsadfs(1,:),[],2);
        traject_bsadf(:,nn_b) = bsadfs'; 
        % bsadfs: vector with the SADFs. This trajectory will be compared with
        % the critical values.
        empirgsadf(nn_b,1) = gsadf;
    if verbose_b == 1;
        toc
        clc
    end
    
    end




end
%parfor nn_b = 1:nrep_b

% Up to here, we have: the boostrapped trayectories, and, the corresponding
% SADF and GSADF for each. 


% Take the quantiles from the bootstrapped distribution, we will use this
% to test the full sample statistics. This will be the first output.
%quantile_sadf = quantile(empirsadf',qe);
%quantile_gsadf = quantile(empirgsadf',qe);
%[criticalvs,~] = obtain_cvs_for_bootstrap(trayect_bsadf,m);

end





function [empirsadf,empirgsadf,true_sadf,true_gsadf] = ...
    bootstrap_from_sample_adjusted_cvs_no_collect(data,nrep_b,T,swindow0,dim,verbose_b,robust) 

% This function computes the bootstrapped distribution for the statistic
% for a given sample data or size T.
% INPUTS
% - data: series (in levels) to bootstrap. Size T
% - nrep_b: number of bootstrap iterations
% - T: sample size of data
% - swindow0: minium number of observations to run the procedure.
% - dim: dim=T-swindow0+1
% - verbose_b = 1 if you want to print time for each iteration

% OUTPUT
% - trayec_bootsrap: a matrix that contains the bootstrapped trayectories
% - empirsadf: the bootstrapped distribution for the SADF
% - empirgsadf: the bootstrapped distribution for the GSADF
% - true_sadf: the value of the SADF in the sample
% - true_gsadf: the value of the GSADF in the sample
%  - collect

% Create space to store some things
%traject_bootstrap = zeros(T,nrep_b); % store bootstraped trajectories
%traject_badf = zeros(dim,nrep_b); % store bootstraped BADF trajectories
%traject_bsadf = zeros(dim,nrep_b); % store bootstraped BSADF trajectories
    if verbose_b ==1
    tic
    end

empirsadf = zeros(nrep_b,1); % store the SADF statistic
empirgsadf = zeros(nrep_b,1); % store the GSADF statistic 

data_init = data(1,1);
ret = data(2:end)-data(1:end-1);
alpha = mean(ret);
resid = ret-alpha;
sigma_hat = std(ret);
r2=(swindow0:1:T)';
% INIT MAIN LOOP
%parfor nn_b = 1:nrep_b
for nn_b = 1:nrep_b

   % Step 1: create bootstraped trayectories as explained in the paper
          if robust ==0
              deltay = alpha+sigma_hat*randn(size(resid));
          else
              deltay = alpha+resid.*randn(size(resid));
          end

   y = [data_init; data_init+cumsum(deltay)] ;  
   %traject_bootstrap(:,nn_b) = y;

   % Compute the statistics
    %THE SUP ADF TEST %%%%%%
    badfs= ADF_FL4(y,swindow0,T);
    sadf=max(badfs);
    %traject_badf(:,nn_b) = badfs;
    empirsadf(nn_b,1)=sadf; 
    % badfs: vector with the ADFs. This trajectory will be compared with
    % the critical values.
%THE GENERALIZED SUP ADF TEST %%%%%%
    bsadfs=zeros(1,dim);
    for v=1:1:size(r2,1);
        rwadft = ADF_FL4_backwards(y(1:r2(v)),swindow0,r2(v));
        bsadfs(1,v)=max(rwadft);
    end;
    gsadf=max(bsadfs(1,:),[],2);
    %traject_bsadf(:,nn_b) = bsadfs'; 
    % bsadfs: vector with the SADFs. This trajectory will be compared with
    % the critical values.
    empirgsadf(nn_b,1) = gsadf;
end

% run it with the true data:

%THE SUP ADF TEST %%%%%%
    true_badfs= ADF_FL4(data,swindow0,T);
    true_sadf=max(true_badfs);
%THE GENERALIZED SUP ADF TEST %%%%%%
    true_bsadfs=zeros(1,dim);
    for v=1:1:size(r2,1)
        rwadft = ADF_FL4_backwards(data(1:r2(v)),swindow0,r2(v));
        true_bsadfs(1,v)=max(rwadft);
    end
    true_gsadf=max(true_bsadfs(1,:),[],2);
    
if verbose_b == 1
    toc
    clc
end

% Up to here, we have: the boostrapped trayectories, and, the corresponding
% SADF and GSADF for each. 


% Take the quantiles from the bootstrapped distribution, we will use this
% to test the full sample statistics. This will be the first output.
%quantile_sadf = quantile(empirsadf',qe);
%quantile_gsadf = quantile(empirgsadf',qe);
%[criticalvs,~] = obtain_cvs_for_bootstrap(trayect_bsadf,m);

end


function [quantile_badfs,quantile_bsadfs,significance_level_badf,significance_level_bsadf] = obtain_cvs_for_bootstrap(trayect_badf,trayect_bsadf,m,target_value)
% THIS FUNCTION OBTAINS THE CRITICAL VALUE TRAJECTORY AND THE ACTUAL
% ADJUSTED SIZE FROM A GIVEN TRAJECTORY OF BSADFs

fun_to_solve_badf = @(x) compute_size(x,trayect_badf,m,target_value);
fun_to_solve_bsadf = @(x) compute_size(x,trayect_bsadf,m,target_value);
actual_size_badf = bisection(fun_to_solve_badf, 0.9, 1);
actual_size_bsadf = bisection(fun_to_solve_bsadf, 0.9, 1);
quantile_badfs = (quantile(trayect_badf',actual_size_badf))';
quantile_bsadfs = (quantile(trayect_bsadf',actual_size_bsadf))';
significance_level_badf = 1-actual_size_badf;
significance_level_bsadf = 1-actual_size_bsadf;
end



function output = compute_size(theo_size,trayect_bsadf,m,target_value)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% - theo_size: the target theoretical size
% - trayect_bsadf a matrix containing BSADFs from bootstrapped previous
% step. Dimensions: dim times nrep, i.e each column is a different
% bootstrapped trajectory.
% - m: minium value of consecutive periods to be considered a bubble. (k in
% the notation of the paper).
% - target_value: target significance value. Usually 0.05
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


series_bsadf = quantile(trayect_bsadf',theo_size)';
bubbles2 = zeros(size(series_bsadf));

for nn_c =1:size(trayect_bsadf,2) %for each trajectory...
bubbles2(2:end,nn_c) = trayect_bsadf(2:end,nn_c)>series_bsadf(2:end);
end


times2 = zeros(size(bubbles2));


for nn=1:size(bubbles2,2)
    for j=2:size(bubbles2,1)    
        if bubbles2(j,nn)== 0
            times2(j,nn)=0;
        else
            times2(j,nn)=times2(j-1,nn)+1;
        end
    
    end
end

% Count the number of detected bubbles per sample defining a bubble as
% having a minimum duration of

det2at12 = (times2==(m*ones(size(times2))));

% here we create a matrix that indicates the number of detected bubbles

numberof2 = sum(det2at12);

% Finally, we count the proportion of samples with at least one of these
% detected bubbles.

FINALcount2 = mean(numberof2>zeros(size(numberof2)));

output = FINALcount2 - target_value;
end



function  [colector_sadf_distribution,colector_gsadf_distribution,colector_sadf_true_values,colector_gsadf_true_values] ...
    = run_double_bootstrap(data_to_use, M,T,swindow0,dim, paralell,verbose_b,robust)

% This function computes the bootstrapped distribution for the statistic
% for a given sample data or size T.
% INPUTS
% - data: series (in levels) to bootstrap. Size T
% - nrep_b: number of bootstrap iterations
% - T: sample size of data
% - swindow0: minium number of observations to run the procedure.
% - dim: dim=T-swindow0+1
% - verbose_b = 1 if you want to print time for each iteration


nrep2 = size(data_to_use,2);
% here we will store the stuff we need:
colector_sadf_distribution = zeros(M,nrep2);
colector_gsadf_distribution = zeros(M,nrep2);
colector_sadf_true_values = zeros(1,nrep2);
colector_gsadf_true_values = zeros(1,nrep2);
%here, each column is the distribution corresponding to a different 
% bootstrap sample

% INIT MAIN LOOP
if paralell == 1
    parfor nn = 1:nrep2
     fprintf('Starting global iter %.0f \n',nn)    
     [empirsadf,empirgsadf,true_sadf,true_gsadf] = ...
    bootstrap_from_sample_adjusted_cvs_no_collect(data_to_use(:,nn),M,T,swindow0,dim,verbose_b,robust); % bootstrap
    % Store results
    colector_sadf_distribution(:,nn) = empirsadf ;
    colector_gsadf_distribution(:,nn) = empirgsadf;
    colector_sadf_true_values(:,nn) = true_sadf ;
    colector_gsadf_true_values(:,nn) = true_gsadf;
    end
else
    for nn = 1:nrep2
     fprintf('Starting global iter %.0f \n',nn)    
     [empirsadf,empirgsadf,true_sadf,true_gsadf] = ...
    bootstrap_from_sample_adjusted_cvs_no_collect(data_to_use(:,nn),M,T,swindow0,dim,verbose_b,robust); % bootstrap
    % Store results
    colector_sadf_distribution(:,nn) = empirsadf ;
    colector_gsadf_distribution(:,nn) = empirgsadf;
    colector_sadf_true_values(:,nn) = true_sadf ;
    colector_gsadf_true_values(:,nn) = true_gsadf;
    end
end
end


function  [cvs_sadf,cvs_gsadf, lambdas_sadf, lambdas_gsadf] = double_bootstrap_analyse_results(sadf_distribution,gsadf_distribution,...
    sadf_true_values,gsadf_true_values,...,
    empirgsadf,empirsadf)
%% Step 1: get lambda such that the true rejection probability is 0.05
true_sizes = [0.1;0.05;0.01];
lambdas_sadf = zeros(size(true_sizes,1),1);
lambdas_gsadf = zeros(size(true_sizes,1),1);

for i = 1:size(true_sizes,1)
target_size = true_sizes(i)
% SADF:
fun_sadf = @(x) compute_true_size2(x,sadf_distribution,sadf_true_values,target_size);
fun_gsadf = @(x) compute_true_size2(x,gsadf_distribution,gsadf_true_values,target_size);
x0 = target_size+0.1;
lambdas_sadf(i) = bisection(fun_sadf,0,x0+0.1);
lambdas_gsadf(i) = bisection(fun_gsadf,0,x0+0.1);
%lambdas_gsadf
end

%% Step 2: use that lambda to get new critical values
% We already loaded a bootrsapped distribution for these. I will load the
% draws I used originally to construct whats on the tables. 

cvs_sadf = quantile(empirsadf,1-lambdas_sadf);
cvs_gsadf = quantile(empirgsadf,1-lambdas_gsadf);

end

function out = compute_true_size2(x_aux,distribution,true_values,target_size)

% x is the nominal level of significance
% distribution is  M (number of bootstrap iters) times nrep (number of
% samples) matrix
% true_values is a row vector that contains the values of the statistic in
% each sample.
% target size: the overall size we are trying to calibrate.
critical_values = quantile(distribution,1-x_aux);
rejection = (true_values>critical_values); 
out = mean(rejection,2)-target_size;
end

function n_bub = compute_number_of_bubbles(traject_bsadf,series_bsadf,m)

bubbles2 = zeros(size(series_bsadf));

for nn_c =1:size(traject_bsadf,2) %for each trajectory...
    bubbles2(2:end,nn_c) = traject_bsadf(2:end,nn_c)>series_bsadf(2:end);
end


times2 = zeros(size(bubbles2));


for nn=1:size(bubbles2,2)
    for j=2:size(bubbles2,1)    
        if bubbles2(j,nn)== 0
            times2(j,nn)=0;
        else
            times2(j,nn)=times2(j-1,nn)+1;
        end
    end
end

% Count the number of detected bubbles per sample defining a bubble as
% having a minimum duration of

n_bub = sum(times2==(m*ones(size(times2))));
end


% function [quantile_badfs,quantile_bsadfs,significance_level_badf,significance_level_bsadf] ...
%     = obtain_cvs_for_bootstrap_precheck(trayect_badf,trayect_bsadf,m,calibrated_cvs_sadf,calibrated_cvs_gsadf,target_value)
% % THIS FUNCTION OBTAINS THE CRITICAL VALUE TRAJECTORY AND THE ACTUAL
% % ADJUSTED SIZE FROM A GIVEN TRAJECTORY OF BSADFs
% 
% 
% 
% fun_to_solve_badf = @(x) compute_size_precheck(x,trayect_badf,calibrated_cvs_sadf,m,target_value);
% fun_to_solve_bsadf = @(x) compute_size_precheck(x,trayect_bsadf,calibrated_cvs_gsadf,m,target_value);
% actual_size_badf = bisection(fun_to_solve_badf, 0.9, 1);
% actual_size_bsadf = bisection(fun_to_solve_bsadf, 0.9, 1);
% quantile_badfs = (quantile(trayect_badf',actual_size_badf))';
% quantile_bsadfs = (quantile(trayect_bsadf',actual_size_bsadf))';
% significance_level_badf = 1-actual_size_badf;
% significance_level_bsadf = 1-actual_size_bsadf;
% end
% 
% function output = compute_size_precheck(theo_size,trayect_bsadf,cvs_full_sample_stat,m,target_value)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Inputs:
% % - theo_size: the target theoretical size
% % - trayect_bsadf a matrix containing BSADFs from bootstrapped previous
% % step. Dimensions: dim times nrep, i.e each column is a different
% % bootstrapped trajectory.
% % - critical value of the full sample statistic used to do the precheck. 
% % basically, if the full value stat is below this critical value, then we
% % don't reject.
% % - m: minium value of consecutive periods to be considered a bubble. (k in
% % the notation of the paper).
% % - target_value: target significance value. Usually 0.05
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% series_bsadf = quantile(trayect_bsadf',theo_size)';
% gsadfs = (max(trayect_bsadf,[],1))';
% bubbles2 = zeros(size(series_bsadf));
% 
% for nn_c =1:size(trayect_bsadf,2) %for each trajectory...
% bubbles2(2:end,nn_c) = trayect_bsadf(2:end,nn_c)>series_bsadf(2:end);
% end
% 
% 
% times2 = zeros(size(bubbles2));
% 
% 
% for nn=1:size(bubbles2,2)
%     for j=2:size(bubbles2,1)    
%         if bubbles2(j,nn)== 0
%             times2(j,nn)=0;
%         else
%             times2(j,nn)=times2(j-1,nn)+1;
%         end
%     end
% end
% 
% % Count the number of detected bubbles per sample defining a bubble as
% % having a minimum duration of
% 
% det2at12 = (times2==(m*ones(size(times2))));
% full_sample = (gsadfs'>cvs_full_sample_stat * ones(1,length(gsadfs)));
% 
% % here we create a matrix that indicates the number of detected bubbles
% 
% numberof2 = sum(det2at12).*full_sample;
% 
% % Finally, we count the proportion of samples with at least one of these
% % detected bubbles.
% 
% FINALcount2 = mean(numberof2>zeros(size(numberof2)));
% 
% output = FINALcount2 - target_value;
% end
%%
function [x,fx,exitFlag] = bisection(f,lb,ub,target,options)
% BISECTION  Fast and robust root-finding method that handles n-dim arrays.
% 
%   [x,fVal,ExitFlag] = BISECTION(f,LB,UB,target,options) finds x within 
%   (LB < x < UB) such that f(x +/- TolX) = target OR f(x) = target +/- TolFun.
% 
%   x = BISECTION(f,LB,UB) finds the root(s) of function f on the interval [LB,
%   UB], i.e. finds x such that f(x) = 0 where LB <= x <= UB. f will never be
%   evaluated outside of the interval specified by LB and UB. f should have only
%   one root and f(UB) and f(LB) must bound it. Elements of x are NaN for
%   instances where a solution could not be found.
% 
%   x = BISECTION(f,LB,UB,target) finds x such that f(x) = target.
% 
%   x = BISECTION(f,LB,UB,target,TolX) will terminate the search when the search
%   interval is smaller than TolX (TolX must be positive).
% 
%   x = BISECTION(f,LB,UB,target,options) solves with the default parameters
%   replaced by values in the structure OPTIONS, an argument created with the
%   OPTIMSET function. Used options are TolX and TolFun. Note that OPTIMSET will
%   not allow arrays for tolerances, so set the fields of the options structure
%   manually for non-scalar TolX or TolFun.
% 
%   [x,fVal] = BISECTION(f,...) returns the value of f evaluated at x.
%
%   [x,fVal,ExitFlag] = BISECTION(...) returns an ExitFlag that describes the
%   exit condition of BISECTION. Possible values of elements of ExitFlag and the
%   corresponding exit conditions are
%
%       1   Search interval smaller than TolX.
%       2   Function value within TolFun of target.
%       3   Search interval smaller than TolX AND function value within 
%           TolFun of target.
%      -1   No solution found.
%      -2   Unbounded root (f(LB) and f(UB) do not span target).
% 
%   Any or all of f(scalar), f(array), LB, UB, target, TolX, or TolFun may be
%   scalar or n-dim arrays. All non-scalar arrays must be the same size. All
%   outputs will be this size.
% 
%   Default values are target = 0, TolX = 1e-6, and TolFun = 0.
% 
%   There is no iteration limit. This is because BISECTION (with a TolX that
%   won't introduce numerical issues) is guaranteed to converge if f is a
%   continuous function on the interval [UB, LB] and f(x)-target changes sign on
%   the interval.
% 
%   The <a href="http://en.wikipedia.org/wiki/Bisection_method">bisection method</a> is a very robust root-finding method. The absolute
%   error is halved at each step so the method converges linearly. However,
%   <a href="http://en.wikipedia.org/wiki/Brent%27s_method">Brent's method</a> (such as implemented in FZERO) can converge
%   superlinearly and is as robust. FZERO also has more features and input
%   checking, so use BISECTION in cases where FZERO would have to be implemented
%   in a loop to solve multiple cases, in which case BISECTION will be much
%   faster because of vectorization.
%
%   Define LB, UB, target, TolX, and TolFun for each specific application using
%   great care for the following reasons:
%     - There is no iteration limit, so given an unsolvable task, BISECTION may
%       remain in an unending loop.
%     - There is no initial check to make sure that f(x) - target changes sign
%       between LB and UB.
%     - Very large or very small numbers can introduce numerical issues.
%
%   Example 1: Find cube root of array 'target' without using NTHROOT and
%   compare speed to using FZERO.
%       options = optimset('TolX', 1e-9);
%       target = [(-100:.1:100)' (-1000:1:1000)'];
% 
%       tic;
%       xfz = zeros(size(target));
%       for ii = 1:numel(target)
%           xfz(ii) = fzero(@(x) x.^3-target(ii), [-20 20], options);
%       end
%       fzero_time = toc
% 
%       tic;
%       xbis = bisection(@(x) x.^3, -20, 20, target, options);
%       bisection_time = toc
% 
%       fprintf('FZERO took %0.0f times longer than BISECTION.\n',...
%                   fzero_time/bisection_time)
% 
%   Example 2: Find roots by varying the function coefficients.
%       [A, B] = meshgrid(linspace(1,2,6), linspace(4,12,10));
%       f = @(x) A.*x.^0.2 + B.*x.^0.87 - 15;
%       xstar = bisection(f,0,5)
% 
%   See also FZERO, FMINBND, OPTIMSET, FUNCTION_HANDLE.
% 
%   [x,fVal,ExitFlag] = BISECTION(f,LB,UB,target,options)

%   FEX URL: http://www.mathworks.com/matlabcentral/fileexchange/28150
%   Copyright 2010-2015 Sky Sartorius
%   Author  - Sky Sartorius
%   Contact - www.mathworks.com/matlabcentral/fileexchange/authors/101715

% Process inputs. 
% Set default values
tolX = 1e-6;
tolFun = 0;
if nargin == 5
    if isstruct(options)
        if isfield(options,'TolX') && ~isempty(options.TolX)
            tolX = options.TolX;
        end 
        if isfield(options,'TolFun') && ~isempty(options.TolFun)
            tolFun = options.TolFun;
        end
    else
        tolX = options;
    end      
end

if (nargin < 4) || isempty(target)
    target = 0;
else
    f = @(x) f(x) - target;
end

ub_in = ub; lb_in = lb; 

% Flip UB and LB if necessary. 
isFlipped = lb > ub;
if any(isFlipped(:))
    ub(isFlipped) = lb_in(isFlipped);
    lb(isFlipped) = ub_in(isFlipped);
    ub_in = ub; lb_in = lb;
end

% Make sure everything is the same size for a non-scalar problem. 
fub = f(ub);
if isscalar(lb) && isscalar(ub)
    % Test if f returns multiple outputs for scalar input.
    if ~isscalar(target)
        id = ones(size(target));
        ub = ub.*id;
        fub = fub.*id;
    elseif ~isscalar(fub)
        ub = ub.*ones(size(fub));
    end
end

% Check if lb and/or ub need to be made into arrays.
if isscalar(lb) && ~isscalar(ub)    
    lb = lb.*ones(size(ub));
elseif ~isscalar(lb) && isscalar(ub)
    id = ones(size(lb));
    ub = ub.*id;
    fub = fub.*id;
end

unboundedRoot = sign(fub).*sign(f(lb)) > 0;
ub(unboundedRoot) = NaN;

% Iterate
ubSign = sign(fub);  
while true
    x = (lb + ub) / 2;
    fx = f(x);
    outsideTolX = abs(ub - x) > tolX;
    outsideTolFun = abs(fx) > tolFun;
    stillNotDone = outsideTolX & outsideTolFun;
    if ~any(stillNotDone(:))
        break;
    end
    select = sign(fx) ~= ubSign;
    lb(select) = x(select);
    ub(~select) = x(~select);
end

% Catch NaN elements of UB, LB, target, or other funky stuff. 
x(isnan(fx)) = NaN;
x(isnan(tolX)) = NaN;
x(isnan(tolFun)) = NaN;
fx(isnan(x)) = NaN;

% Characterize results. 
if nargout > 1 && nnz(target(:))
    fx = fx + target;
end
if nargout > 2 
    exitFlag                                    = +~outsideTolX;
    exitFlag(~outsideTolFun)                    =  2;
    exitFlag(~outsideTolFun & ~outsideTolX)     =  3;
    exitFlag(isnan(x))                          = -1;
    exitFlag(unboundedRoot)                     = -2;
end

end
