%Function to apply GLM-PAC procedure.
%
%INPUTS.
%  Vlo = low frequency band signal.
%  Vhi = high frequency band signal.
%  nCtlPts = the number of control points to use in the spline fitting of phase.
%
%OPTIONAL INPUT.
%  The 4th input is optional.  Set the 4th input to:
%    'noplot' to prevent plotting of results.
%    'AIC'    to compute the # control points using AIC.  When this option is used,
%             the number of control points that minimizes the AIC is used, and the
%             3rd input (nCtlPts) is ignored.
%
%OUTPUTS.
%  PAC   = the "r-value" for the PAC-AAC model.
%  PAC_CI = the 95% confidence intervals for PAC.
%  AAC    = estimate of the effect of ampLO on ampHi.
%  AAC_CI = the 95% confidence intervals for AAC.
%  AAC_p  = p-value for AAC.
%  nCtlPts = the number of control points used.
%
%  By default, this function plots the results.

function [PAC,PAC_CI, AAC,AAC_CI,AAC_p, nCtlPts] = GLM_PAC_AAC(Vlo, Vhi, nCtlPts, varargin)

  %Compute phase and amplitude.
  phi = angle(hilbert(Vlo));
  amp = abs(hilbert(Vhi));
  ampLO = abs(hilbert(Vlo));

  %Compute AIC to determine number of control points.
  
  if ~isempty(varargin) & strcmp(varargin{1}, 'AIC')
      fprintf('Running AIC ... \n')
      Y = amp';
      CtlPts = (4:1:30);
      %Compute the AIC.
      AIC = zeros(size(CtlPts));
      for k=1:length(CtlPts)
          nCtlPts = CtlPts(k);
          X1 = spline_phase0(phi',nCtlPts);
          
          [b1, dev, stats1]    = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off');
          
          AIC(k) = dev + 2*nCtlPts;
          fprintf([num2str(nCtlPts) ' ' num2str(AIC(k)) '\n'])
      end
      %Select the # control points from AIC, and plot the AIC.
      [~, imn] = min(AIC);
      nCtlPts = CtlPts(imn);
      figure(1);  clf();
      plot(CtlPts, AIC, 'k', 'LineWidth', 2)
      hold on
      plot([nCtlPts, nCtlPts], [min(AIC) max(AIC)], 'r', 'LineWidth', 2)
      hold off
      axis tight
      xlabel('# control points')
      ylabel('AIC')
      fprintf(['Suggested number of control points is ' num2str(nCtlPts) '\n'])
  end
  
  %Define variables for GLM procedure.
  Y = amp'; 
  X1 = spline_phase0(phi',nCtlPts);
  X2 = [X1,ampLO',ampLO'.*sin(phi'),ampLO'.*cos(phi')]; % the model takes into account the low-frequency amplitude
  XC = ones(size(Y));  

  %Perform GLM.
  [b1,  ~, stats1]  = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off');
  [b2,  ~, stats2]  = glmfit(X2,  Y, 'gamma', 'link', 'log', 'constant', 'off');
  [bC,  ~, statsC]  = glmfit(XC,  Y, 'gamma', 'link', 'log', 'constant', 'off');

  %Define dense phase points for interpolation.
  phi0 = linspace(-pi,pi,100);
  X0 = spline_phase0(phi0',nCtlPts);
  X2eval = [X0,exp(bC)*ones(size(phi0))', ...
               exp(bC).*sin(phi0)', ...
               exp(bC).*cos(phi0)']; %evaluate on some fixed amplitude

  %Determine spline fit and CI.
  [spline0, dylo, dyhi] = glmval(b1,X0,'log',stats1,'constant', 'off');
  splineU = spline0+dyhi;
  splineL = spline0-dylo;
  
  [spline2, dylo2, dyhi2] = glmval(b2,X2eval,'log',stats2,'constant', 'off');
  splineU2 = spline2+dyhi2;
  splineL2 = spline2-dylo2;

  %Determine null fit and CI.
  [null0, dylo, dyhi] = glmval(bC,ones(size(phi0)),'log',statsC,'constant', 'off');
  nullU = null0+dyhi;
  nullL = null0-dylo;

  %Find the max absolute percentage change between the two models.
  [r,imx]   = max(abs(1-spline0./null0));
  [r2,imx2] = max(abs(1-spline2./null0));

  %Determine CI for the measure r.
  M = 10000;
  bMC = b1*ones(1,M) + sqrtm(stats1.covb)*normrnd(0,1,nCtlPts,M);
  splineMC = glmval(bMC,X0,'log',stats1,'constant', 'off');
  nullMC   = exp(bC)*ones(1,M);%mean(splineMC,1);
  mx = zeros(M,1);
  for k=1:M
      mx(k) = max(abs(1-splineMC(:,k)./nullMC(k)));
  end
  r_CI = quantile(mx, [0.025, 0.975]);
  
  %and for r2
  M = 10000;
  bMC = b2*ones(1,M) + sqrtm(stats2.covb)*normrnd(0,1,nCtlPts+3,M);
%   splineMC = zeros(length(phi0),M);
%   for j=1:M
%       bMC = b2 + sqrtm(stats2.covb)*normrnd(0,1,nCtlPts+1,1);
%       X2eval = [X0,randsample(ampLO,1)*ones(size(phi0))'];
%       splineMC(:,j) = glmval(bMC,X2eval,'log',stats2,'constant', 'off');
%   end
  splineMC = glmval(bMC,X2eval,'log',stats2,'constant', 'off');
  nullMC   = exp(bC)*ones(1,M); %mean(splineMC,1);
  mx = zeros(M,1);
  for k=1:M
      mx(k) = max(abs(1-splineMC(:,k)./nullMC(k)));
  end
  r2_CI = quantile(mx, [0.025, 0.975]);

  %Plot the results.
  if isempty(varargin) || ~strcmp(varargin{1}, 'noplot')
      %figure(2)
%       plot(phi0, spline0, 'r', 'LineWidth', 2)
%       hold on
%       plot(phi0, splineU, ':r', 'LineWidth', 2)
%       plot(phi0, splineL, ':r', 'LineWidth', 2)
      
      plot(phi0, null0, 'k', 'LineWidth', 2);  hold on
      plot(phi0, nullL, 'k:', 'LineWidth', 2)
      plot(phi0, nullU, 'k:', 'LineWidth', 2)
      
      plot(phi0,spline2,'b','LineWidth',2)
      plot(phi0,splineU2,'b:','LineWidth',2)
      plot(phi0,splineL2,'b:','LineWidth',2)
      
      plot([phi0(imx2) phi0(imx2)], [null0(imx2), spline2(imx2)], 'LineWidth', 2)
      %plot([phi0(imx) phi0(imx)], [null0(imx), spline0(imx)], 'LineWidth', 2)
      hold off
      ylabel('Amplitude')
      xlabel('Phase')
      axis tight
  end
  
  PAC    = r2;
  PAC_CI = r2_CI;
  AAC    = exp(b2(nCtlPts+1));
  AAC_CI = [exp(b2(nCtlPts+1)-2*stats2.se(nCtlPts+1)), ...
            exp(b2(nCtlPts+1)+2*stats2.se(nCtlPts+1))];
  AAC_p  = stats2.p(nCtlPts+1);
  
%   fprintf(['PAC: r = ' num2str(r2,3) ' [' num2str(r2_CI(1),3) '-' num2str(r2_CI(2),3) '] \n'])
%   
%   fprintf(['AAC: modulation of A_high by A_lo * ' num2str(AAC,3) ...
%             ' [' num2str(AAC_CI(1),3) '-' ...
%                  num2str(AAC_CI(2),3) ']' ...
%                 ', p=' num2str(AAC_p,3) ...
%             '\n'])
  

        
        
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
