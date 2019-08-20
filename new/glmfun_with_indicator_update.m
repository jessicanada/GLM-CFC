function [XX,P,I] = glmfun_with_indicator_update(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,pval,ci,varargin)
% a modification of glmfun, allowing the user to include an indicator
% function to the model, differentiating between a 'pre' and 'post'
% condition (e.g. before and after stimulus, with or without drug), to test
% for the effect of condition on PAC
%
% INPUTS:
% Vlo_pre:                    Low frequency signal from 'pre' condition
% Vhi_pre:                    High frequency signal from 'pre' condition
% Vlo_post:                   Low frequency signal from 'post' condition
% Vhi_post:                   High frequency signal from 'post' condition
% pval:                       'theoretical' gives analytic p-values for R
%                             'empirical' gives bootstrapped p-values for R
% ci:                         'ci' gives confidence intervals, 'none' gives no confidence intervals
% varargin:                   optionally, include the parameter q indicating which quantiles
%                             of AmpLo you'd like to fit over
%
% OUTPUTS:
% XX.rpac:      R_PAC value, confidence intervals XX.rPAC_CI
% XX.raac:      R_AAC value, confidence intervals XX.rAAC_CI
% XX.rcfc:      R_CFC value, confidence intervals XX.rCFC_CI
% XX.null:      3D surface for null model in Phi_low, A_low, A_high space
% XX.PAC:       3D surface for PAC model in Phi_low, A_low, A_high space
% XX.AAC:       3D surface for AAC model in Phi_low, A_low, A_high space
% XX.CFC:       3D surface for CFC model in Phi_low, A_low, A_high space
% P.rpac:       p-value for RPAC statistic
% P.raac:       p-value for RAAC statistic
% P.rcfc:       p-value for RCFC statistic
% I.pval:       p-value for indicator variable: low value implies significant effect of
%               indicator on PAC
% I.beta:       output beta values from full model including indicator
%               variable
      
  nCtlPts = 10;
  
  %Compute phase and amplitude.
  phi = [angle(hilbert(Vlo_pre));angle(hilbert(Vlo_post))];
  amp = [abs(hilbert(Vhi_pre));abs(hilbert(Vhi_post))];
  ampLO = [abs(hilbert(Vlo_pre)); abs(hilbert(Vlo_post))];
  POST = [zeros(size(Vlo_pre)); ones(size(Vlo_post))]; % (0) 1 indicates (not) in Post scenario 
  
  %Define variables for GLM procedure.
  Y = amp;                                                                 %high frequency amplitude
  X1 = spline_phase0(phi',nCtlPts);                                         %low frequency phase
  X2 = [ones(size(Y)),ampLO];                                              %low frequency amplitude
  X3 = [X1,ampLO,sin(phi).*ampLO,cos(phi).*ampLO];                     %low frequency phase, low frequency amplitude, interaction terms
  X4 = [X3, X1.*POST, ampLO.*POST];
  X5 = [ones(size(Y)),ampLO,ampLO.*POST,POST];
  X6 = [X1, X1.*POST];
  XC = ones(size(Y));                                                       %ones (null)

  %Perform GLM.
  [b1,  dev1, stats1]  = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off'); %Phi_low
  [b2,  dev2, stats2] = glmfit(X2, Y, 'gamma','link','log','constant','off');       %A_low
  [b3,  dev3, stats3] = glmfit(X3, Y, 'gamma','link','log','constant','off');       %Phi_low, A_low
  [bC, dev0, statsC] = glmfit(XC, Y, 'gamma', 'link', 'log', 'constant', 'off');    %null
  [b4, dev4, stats4] = glmfit(X4, Y, 'gamma', 'link', 'log', 'constant', 'off');    %Phi_low, A_low, condition
  [b5, dev5, stats5] = glmfit(X5, Y, 'gamma', 'link', 'log', 'constant', 'off');    %A_low, condition
  [b6, dev6, stats6] = glmfit(X6, Y, 'gamma', 'link','log','constant','off');       %Phi_low, condition
  
  %Chi^2 test between nested models (theoretical p-values)
  chi0 = 1-chi2cdf(dev0-dev3,12);
  chi1 = 1-chi2cdf(dev1-dev3,3); %Between PAC and PACAAC Model, if low AAC is present
  chi2 = 1-chi2cdf(dev2-dev3,11); %Between AAC and PACAAC Model, if low PAC is present
  
  X7 = [X3,X1.*POST]; %model 4 without additional A_low term
  X8 = [X3,ampLO.*POST]; %model 4 without additional Phi_low term
  [~,dev7,~] = glmfit(X7,Y,'gamma','link','log','constant','off');
  [b8,dev8,stats8] = glmfit(X8,Y,'gamma','link','log','constant','off');
  
  I.pval_AAC = 1-chi2cdf(dev7-dev4,1);
  I.pval_PAC = 1-chi2cdf(dev8-dev4,10);
  [~,pvalue] = lratiotest(dev8,dev4,10);
  I.pval_PAC_lrtest = pvalue;
  I.beta = b4;
  
 %create 3d model surfaces
      phi0 = linspace(-pi,pi,100);
      ampSORT = sort(ampLO);
      if ~isempty(varargin)
          q = varargin{1};
          ind = find(ampLO>quantile(ampLO,q) & ampLO<quantile(ampLO,1-q));
          ampSORT = sort(ampLO(ind));
      end
      XXC = []; XX1 = []; XX2 = []; XX3 = []; XX4 = []; XX5 = [];XX8=[];
      ampAXIS = [];
      YC = ones(size(phi0));                                            %null model
      Y1 = spline_phase0(phi0',nCtlPts);                                %PAC model, function of phiLo
      [splineC, ~, ~] = glmval(bC,YC,'log',statsC,'constant', 'off');   %null
      [spline1, ~, ~] = glmval(b1,Y1,'log',stats1,'constant', 'off');   %PAC
      count = 1;
      for i = 1:100:length(ampSORT)                                       %fit model over values of Alow
          Y3 = [Y1,ampSORT(i)*ones(size(phi0))',ampSORT(i)*sin(phi0'),ampSORT(i)*cos(phi0')];        %CFC model, function of phiLo, ampLo, phiLo*ampLo
          [spline3, ~, ~] = glmval(b3,Y3,'log',stats3,'constant', 'off');
          XX3(:,count) = spline3;
          for j = 1:2
              POST_val = ones(size(Y1))*(j-1); %all 1s or all 0s
              Y4 = [Y3,Y1.*POST_val,ampSORT(i).*ones(100,1)*(j-1)];
              [spline4,~,~] = glmval(b4,Y4,'log',stats4,'constant','off');
              XX4(:,count,j) = spline4; %4D
              Y8 = [Y3,ampSORT(i).*ones(100,1)*(j-1)];
              [spline8,~,~] = glmval(b8,Y8,'log',stats8,'constant','off');
              XX8(:,count,j) = spline8;
          end
          ampAXIS(count) = ampSORT(i);
          count = count+1;
      end
      L = length(1:100:length(ampSORT));
      XX1 = repmat(spline1,1,L); %PAC model constant in PhiLow dimension
      XXC = repmat(splineC,1,L); %null model constant in PhiLow, Alow dimensions
      temp = ampSORT(1:100:end);
      Y2 = [ones(size(temp)),temp];
      [spline2,~,~] = glmval(b2,Y2,'log',stats2,'constant','off'); %AAC model, function of Alow
      Xtemp = repmat(spline2,1,100);                                         %AAC model constant in PhiLow dimension
      for l = 1:2
        POST_val = ones(size(temp))*(l-1);
        Y5 = [ones(size(temp)),temp,POST_val.*temp,POST_val];
        [spline5,~,~] = glmval(b5,Y5,'log',stats5,'constant','off');
        XX5(:,:,l) = repmat(spline5,1,100)';
      end
      
      XX.AAC = Xtemp'; XX.null = XXC; XX.PAC = XX1;XX.CFC = XX3; XX.CFC_indicator = XX4;         %3D model surfaces
      XX.AAC_indicator = XX5; XX.model4 = XX4; XX.model8 = XX8;
      
      XX.ampAXIS = ampAXIS; XX.phi0 = phi0;     %axes
      XX.rpac = max(max((abs(1-XX1./XXC))));
      XX.rpac_new = max(max(abs(1-Xtemp'./XX3)));
      XX.raac = max(max((abs(1-Xtemp'./XXC))));
      XX.raac_new = max(max((abs(1-XX1./XX3))));
      XX.rcfc = max(max((abs(1-XX3./XXC))));
      
      % calculate new R values
      
      AAC_4D = repmat(Xtemp',1,1,2);
      
      XX.AAC_4D = AAC_4D;
      XX.rpac_new_indicator = max(max(max(abs(1-XX5./XX4))));
      

  
  if exist('pval','var') && strcmp(pval, 'empirical')
    M = minvals(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post); %find empirical p-values
    P.rpac = max(.5,length(find(M.rpac>XX.rpac)))/length(M.rpac);
    P.raac = max(.5,length(find(M.raac>XX.raac)))/length(M.raac);
    P.rcfc = max(.5,length(find(M.rcfc>XX.rcfc)))/length(M.rcfc);
    P.rpac_indicator = max(.5,length(find(M.rpac_indicator>XX.rpac_new_indicator)))/length(M.rpac_indicator);
  elseif exist('pval','var') && strcmp(pval, 'theoretical')
    P.rpac = chi2;  %use theoretical p-values
    P.raac = chi1;
    P.rcfc = chi0;
  else
    P = 'No p-values output';
  end
  
  %distance distributions between model 4 and model 8
  M = 1000;
  phi0 = linspace(-pi,pi,100);                                    %null model
  Y1 = spline_phase0(phi0',nCtlPts);    
         
  bMC_w = b4*ones(1,M) + sqrtm(stats4.covb)*normrnd(0,1,length(b4),M);
  bMC_wo = b8*ones(1,M) + sqrtm(stats8.covb)*normrnd(0,1,length(b8),M);
  
  d = zeros(1,M);
  for k = 1:M
      XX_w = []; XX_wo = [];
      b_w = bMC_w(:,k); b_wo = bMC_wo(:,k);
      count = 1;
      for i = 1:100:length(ampSORT)                                       %fit model over values of Alow
          Y3 = [Y1,ampSORT(i)*ones(size(phi0))',ampSORT(i)*sin(phi0'),ampSORT(i)*cos(phi0')];        %CFC model, function of phiLo, ampLo, phiLo*ampLo
          for j = 1:2
              POST_val = ones(size(Y1))*(j-1); %all 1s or all 0s
              Y4 = [Y3,Y1.*POST_val,ampSORT(i).*ones(100,1)*(j-1)];
              [spline_w,~,~] = glmval(b_w,Y4,'log',stats4,'constant','off');
              XX_w(:,count,j) = spline_w;
              %Y8 = [Y3,ampSORT(i).*ones(100,1)*(j-1)];
              %[spline_wo,~,~] = glmval(b_wo,Y8,'log',stats8,'constant','off');
              %XX_wo(:,count,j) = spline_wo;
          end
          count = count+1;
      end
      d(k) = max(max(max((abs(1-XX_w./XX4)))));
  end

  d_4_8 = max(max(max((abs(1-XX8./XX4)))));
  I.d_pval = max(.5,length(find(d>d_4_8)))/length(d);
  
  %confidence intervals
  if exist('ci','var') && strcmp(ci, 'ci')
    
    phi0 = linspace(-pi,pi,100);
    X0 = spline_phase0(phi0',nCtlPts);
    Amax = max(ampSORT); Amin = min(ampSORT); stepsize = (Amax-Amin)/99;
    X2eval = Amin:stepsize:Amax; %evaluate on all amplitudes
    X2eval = [ones(size(phi0))',X2eval'];
    
    M = 10000;
    
    %R_PAC ci
    bMC = b2*ones(1,M) + sqrtm(stats2.covb)*normrnd(0,1,2,M);
    splineAAC = glmval(bMC,X2eval,'log',stats2,'constant', 'off');
     
    [m] = max(abs(1-Xtemp'./XX3)); %find point of maximum distance between A_low and A_low,Phi_low models
    [~,j] = max(m);                         %j ampLO index
    [m2] = max(abs(1-Xtemp'./XX3),[],2);
    [~,j2] = max(m2);                       %j2 phiLO index
    phiLOW = phi0(j2);
    
    bMC = b3*ones(1,M) + sqrtm(stats3.covb)*normrnd(0,1,nCtlPts+3,M);
    %A_low,Phi_low model at fixed Phi_low
    Y_phi = spline_phase0(phiLOW*ones(size(phi0)),nCtlPts); %model 1, function of phiLo
    AmpLOW = [Amin:stepsize:Amax]';
    Y_phi = [Y_phi,AmpLOW]; %model 2, function of phiLo, ampLo
    Y_phi = [Y_phi,AmpLOW*sin(phiLOW'),AmpLOW*cos(phiLOW')];
    spline_phi = glmval(bMC,Y_phi,'log',stats3,'constant', 'off'); %CFC model, constant in phi_low dim
    %A_low,Phi_low model at fixed A_low
    Y_a = spline_phase0(phi0',nCtlPts); %model 1, function of phiLo
    Y_a = [Y_a,ampAXIS(j)*ones(size(phi0))']; %model 2, function of phiLo, ampLo
    Y_a = [Y_a,ampAXIS(j)*sin(phi0'),ampAXIS(j)*cos(phi0')];
    spline_a = glmval(bMC,Y_a,'log',stats3,'constant','off'); %CFC model, constant in A_low dim
    
    %A_low model at fixed A_low
    bMC = b2*ones(1,M) + sqrtm(stats2.covb)*normrnd(0,1,2,M);
    X2eval_fixed_A = [ones(size(phi0))',ones(size(phi0))'*ampAXIS(j)];
    splineMC = glmval(bMC,X2eval_fixed_A,'log',stats2,'constant', 'off');
    splineAAC_a = splineMC;
      
    mx = zeros(M,1);
      for k = 1:M
          m_a = max(abs(1-splineAAC_a(:,k)./spline_a(:,k))); %fixed A
          m_phi = max(abs(1-splineAAC(:,k)./spline_phi(:,k))); %fixed Phi
          mx(k) = max(m_a,m_phi);
      end
      XX.rpac_ci = quantile(mx,[0.025,0.975]);
  
  
    %R_AAC ci
    bMC = b1*ones(1,M) + sqrtm(stats1.covb)*normrnd(0,1,nCtlPts,M);
    splinePAC = glmval(bMC,X0,'log',stats1,'constant', 'off');
    
    [m] = max(abs(1-XX1./XX3)); %find point of maximum distance between A_low and A_low,Phi_low models
    [~,j] = max(m);                         %j ampLO index
    [m2] = max(abs(1-XX1./XX3),[],2);
    [~,j2] = max(m2);                       %j2 phiLO index
    phiLOW = phi0(j2);
    
    bMC = b3*ones(1,M) + sqrtm(stats3.covb)*normrnd(0,1,nCtlPts+3,M);
    %A_low,Phi_low model at fixed Phi_low
    Y_phi = spline_phase0(phiLOW*ones(size(phi0)),nCtlPts);
    AmpLOW = [Amin:stepsize:Amax]';
    Y_phi = [Y_phi,AmpLOW];
    Y_phi = [Y_phi,AmpLOW*sin(phiLOW'),AmpLOW*cos(phiLOW')];
    spline_phi = glmval(bMC,Y_phi,'log',stats3,'constant', 'off'); %A_low,Phi_low model, constant in phi_low dim
    %A_low,Phi_low model at fixed A_low
    Y_a = spline_phase0(phi0',nCtlPts); 
    Y_a = [Y_a,ampAXIS(j)*ones(size(phi0))']; 
    Y_a = [Y_a,ampAXIS(j)*sin(phi0'),ampAXIS(j)*cos(phi0')];
    spline_a = glmval(bMC,Y_a,'log',stats3,'constant','off');
    %Phi_low model at fixed Phi_low
    X0_fixed_phi = spline_phase0(phiLOW*ones(size(phi0)),nCtlPts);
    bMC = b1*ones(1,M) + sqrtm(stats1.covb)*normrnd(0,1,nCtlPts,M);
    splineMC = glmval(bMC,X0_fixed_phi,'log',stats1,'constant', 'off');
    splinePAC_phi = splineMC;
      
    mx = zeros(M,1);
      for k = 1:M
          m_a = max(abs(1-splinePAC(:,k)./spline_a(:,k))); %fixed A
          m_phi = max(abs(1-splinePAC_phi(:,k)./spline_phi(:,k))); %fixed Phi
          mx(k) = max(m_a,m_phi);
      end
      XX.raac_ci = quantile(mx,[0.025,0.975]);
  end

end

% Bootstrapped p-values
function M = minvals(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post)
    K = 100;
    RPAC = zeros(1,K); RCFC = zeros(1,K); RAAC = zeros(1,K);  RPAC_indicator = zeros(1,K);
    for i = 1:K
        i
        Vhi = [Vhi_pre;Vhi_post];
        Vhi_prime = AAFT(Vhi,1);
        Vhi_prime_pre = Vhi_prime(1:length(Vhi_pre));
        Vhi_prime_post = Vhi_prime(length(Vhi_pre)+1:end);
        [XX] = glmfun_with_indicator(Vlo_pre,Vlo_post,Vhi_prime_pre,Vhi_prime_post,'none','none',.05);
        RPAC(i) = XX.rpac;
        RCFC(i) = XX.rcfc;
        RAAC(i) = XX.raac;
        RPAC_indicator(i) = XX.rpac_new_indicator;
    end
    M.rpac = RPAC; M.rcfc = RCFC; M.raac = RAAC; M.rpac_indicator = RPAC_indicator;
end

% Generate a design matrix X (n by nCtlPts) for a phase signal (n by 1)
function X = spline_phase0(phase,nCtlPts)

% Define Control Point Locations
c_pt_times_all = linspace(0,2*pi,nCtlPts+1);

s = 0.5;  % Define Tension Parameter

% Construct spline regressors
X = zeros(length(phase),nCtlPts);
for i=1:length(phase)  
    nearest_c_pt_index = max(find(c_pt_times_all<=mod(phase(i),2*pi)));
    nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
    next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
    u = (mod(phase(i),2*pi)-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    X(i,mod(nearest_c_pt_index-2:nearest_c_pt_index+1,nCtlPts)+1) = p;   
end

end

