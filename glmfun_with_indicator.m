function [XX,P,I] = glmfun_with_indicator(Vlo_pre,Vlo_post,Vhi_pre,Vhi_post,pval,ci,varargin)
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
  POST = [zeros(size(Vlo_pre)); ones(size(Vlo_post))];
  
  %Define variables for GLM procedure.
  Y = amp;                                                                 %high frequency amplitude
  X1 = spline_phase0(phi',nCtlPts);                                         %low frequency phase
  X2 = [ones(size(Y)),ampLO];                                              %low frequency amplitude
  X3 = [X1,ampLO,sin(phi).*ampLO,cos(phi).*ampLO];                     %low frequency phase, low frequency amplitude, interaction terms
  X4 = [X3,X1.*POST];
  XC = ones(size(Y));                                                       %ones (null)

  %Perform GLM.
  [b1,  dev1, stats1]  = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off'); %PAC
  [b2,  dev2, stats2] = glmfit(X2, Y, 'gamma','link','log','constant','off');       %AAC
  [b3,  dev3, stats3] = glmfit(X3, Y, 'gamma','link','log','constant','off');       %CFC
  [bC, dev0, statsC] = glmfit(XC, Y, 'gamma', 'link', 'log', 'constant', 'off');    %null
  [b4, dev4, stats4] = glmfit(X4, Y, 'gamma', 'link', 'log', 'constant', 'off');
  
  %Chi^2 test between nested models (theoretical p-values)
  chi0 = 1-chi2cdf(dev0-dev3,12);
  chi1 = 1-chi2cdf(dev1-dev3,3); %Between PAC and PACAAC Model, if low AAC is present
  chi2 = 1-chi2cdf(dev2-dev3,11); %Between AAC and PACAAC Model, if low PAC is present
  chi3 = 1-chi2cdf(dev3-dev4,10);
  
  I.pval = chi3;
  I.beta = b4;
  
 %create 3d model surfaces
      phi0 = linspace(-pi,pi,100);
      ampSORT = sort(ampLO);
      if ~isempty(varargin)
          q = varargin{1};
          ind = find(ampLO>quantile(ampLO,q) & ampLO<quantile(ampLO,1-q));
          ampSORT = sort(ampLO(ind));
      end
      XXC = []; XX1 = []; XX2 = []; XX3 = [];
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
      
      XX.AAC = Xtemp'; XX.null = XXC; XX.PAC = XX1;XX.CFC = XX3;          %3D model surfaces
         
      XX.ampAXIS = ampAXIS; XX.phi0 = phi0;     %axes
      XX.rpac = max(max((abs(1-XX1./XXC))));
      XX.raac = max(max((abs(1-Xtemp'./XXC))));
      XX.rcfc = max(max((abs(1-XX3./XXC))));

  
  if exist('pval','var') && strcmp(pval, 'empirical')
    M = minvals(Vlo,Vhi); %find empirical p-values
    P.rpac = max(.5,length(find(M.rpac>XX.rpac)))/length(M.rpac);
    P.raac = max(.5,length(find(M.raac>XX.raac)))/length(M.raac);
    P.rcfc = max(.5,length(find(M.rcfc>XX.rcfc)))/length(M.rcfc);
  elseif exist('pval','var') && strcmp(pval, 'theoretical')
    P.rpac = chi2;  %use theoretical p-values
    P.raac = chi1;
    P.rcfc = chi0;
  else
    P = 'No p-values output';
  end
  
  if exist('ci','var') && strcmp(ci, 'ci')
        phi0 = linspace(-pi,pi,100);
        X0 = spline_phase0(phi0',nCtlPts);
        Amax = max(ampSORT); Amin = min(ampSORT); stepsize = (Amax-Amin)/99;
        X2eval = Amin:stepsize:Amax; %evaluate on all amplitudes
        X2eval = [ones(size(phi0))',X2eval'];
      %Determine CI for the measure RPAC.
      M = 10000;
      bMC = b1*ones(1,M) + sqrtm(stats1.covb)*normrnd(0,1,nCtlPts,M);
      splineMC = glmval(bMC,X0,'log',stats1,'constant', 'off');
      mx = zeros(M,1);
      for k=1:M
          mx(k) = max(abs(1-splineMC(:,k)./splineC));
      end
      r_CI = quantile(mx, [0.025, 0.975]);
      XX.rpac_ci = r_CI;

      %and for rAAC
      M = 10000;
      bMC = b2*ones(1,M) + sqrtm(stats2.covb)*normrnd(0,1,2,M);
      splineMC = glmval(bMC,X2eval,'log',stats2,'constant', 'off');
      mx = zeros(M,1);
      for k=1:M
          mx(k) = max(abs(1-splineMC(:,k)./splineC));
      end
      r2_CI = quantile(mx, [0.025, 0.975]);
      XX.raac_ci = r2_CI;

      %and for rCFC
      M = 10000;
      [m,~] = max(abs(1-XX3./XXC)); %find point of maximum distance between null, CFC models
      [~,j] = max(m);                         %j ampLO, I(j) phiLO
      bMC = b3*ones(1,M) + sqrtm(stats3.covb)*normrnd(0,1,nCtlPts+3,M);
      Y1 = spline_phase0(phi0',nCtlPts); %model 1, function of phiLo
      Y2 = [Y1,ampAXIS(j)*ones(size(phi0))']; %model 2, function of phiLo, ampLo
      Y3 = [Y2,ampAXIS(j)*sin(phi0'),ampAXIS(j)*cos(phi0')];
      splineMC = glmval(bMC,Y3,'log',stats3,'constant', 'off');
      mx = zeros(M,1);
      for k=1:M
          mx(k) = max(abs(1-splineMC(:,k)./splineC));
      end
      r3_CI = quantile(mx, [0.025, 0.975]);
      XX.rcfc_ci = r3_CI;
  end

end

% Bootstrapped p-values
function M = minvals(Vlo,Vhi)
    K = 100;
    RPAC = zeros(1,K); RCFC = zeros(1,K); RAAC = zeros(1,K);
    N = zeros(1,K); L = zeros(1,K);
    for i = 1:K
        Vhi_prime = AAFT(Vhi,1);
        [XX] = glmfun(Vlo,Vhi_prime','none');        %compute R statistics between Vhi and shifted Vlo
        RPAC(i) = XX.rpac;
        RCFC(i) = XX.rcfc;
        RAAC(i) = XX.raac;
    end
    M.rpac = RPAC; M.rcfc = RCFC; M.raac = RAAC; M.shiftN = N; M.shiftL = L;
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

