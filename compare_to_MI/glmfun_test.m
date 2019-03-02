function [R,CHI,XX,bC] = glmfun_test(Vlo, Vhi, nCtlPts,varargin)
%INPUTS:
%Vlo: Low frequency signal (4-7Hz)
%Vhi: High frequency signal (100-140Hz)
%nCtlPts: number of control points, for spline phase

%OUTPUTS:
%XX: Structure of spline matrices, each column corresponding to a different
%ampLo value
%ri: maximum difference from model i to the null model
%ri_CI: confidence intervals for the value ri
%chitest: Chi^2 test, determining the significance of the interaction term
%ampLO*phiLO in model 3
      
  %Compute phase and amplitude.
  phi = angle(hilbert(Vlo));
  amp = abs(hilbert(Vhi));
  ampLO = abs(hilbert(Vlo));
  AFIX = max(ampLO);
  
  %Define variables for GLM procedure.
  Y = amp'; 
  X1 = spline_phase0(phi',nCtlPts);
  X2 = [ones(size(Y)),ampLO']; % the model takes into account the low-frequency amplitude
  X3 = [X1,ampLO',sin(phi').*ampLO',cos(phi').*ampLO'];
  XC = ones(size(Y));  

  %Perform GLM.
  [b1,  dev1, stats1]  = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off');
  [b2,  dev2, stats2] = glmfit(X2, Y, 'gamma','link','log','constant','off');
  [b3,  dev3, stats3] = glmfit(X3, Y, 'gamma','link','log','constant','off');
  [bC, dev0, statsC] = glmfit(XC, Y, 'gamma', 'link', 'log', 'constant', 'off');
  
  
  chi0 = 1-chi2cdf(dev0-dev3,12);
  chi1 = 1-chi2cdf(dev1-dev3,3); %Between PAC and PACAAC Model, if low AAC is present
  chi2 = 1-chi2cdf(dev2-dev3,11); %Between AAC and PACAAC Model, if low PAC is present
  CHI.cfcnull = chi0;
  CHI.cfcpac = chi1;
  CHI.cfcaac = chi2;
  
 
      phi0 = linspace(-pi,pi,100);
      ampSORT = sort(ampLO);
      XXC = []; XX1 = []; XX2 = []; XX3 = [];
      ampAXIS = [];
      YC = ones(size(phi0)); %null model
      Y1 = spline_phase0(phi0',nCtlPts); %model 1, function of phiLo
      [splineC, ~, ~] = glmval(bC,YC,'log',statsC,'constant', 'off');
      [spline1, ~, ~] = glmval(b1,Y1,'log',stats1,'constant', 'off');
      count = 1;
      for i = 1:100:length(ampLO) %CHANGED 10 to 100
          ampSORT(i) = exp(bC);
          Y2 = [Y1,ampSORT(i)*ones(size(phi0))']; %model 2, function of phiLo, ampLo
          Y3 = [Y2,ampSORT(i)*sin(phi0'),ampSORT(i)*cos(phi0')]; %model 4, function of phiLo, ampLo, phiLo*ampLo
          [spline3, ~, ~] = glmval(b3,Y3,'log',stats3,'constant', 'off');
          XX3(:,count) = spline3;
          ampAXIS(count) = ampSORT(i);
          count = count+1;
      end
      L = length(1:100:length(ampLO));
      XX1 = repmat(spline1,1,L); %CHANGED 6400 to 640
      XX1 = XX1;
      XXC = repmat(splineC,1,L); %CHANGED 6400 to 640
      XXC = XXC;
      temp = ampSORT(1:100:end);
      Y2 = [ones(size(temp')),temp'];
      [spline2,~,~] = glmval(b2,Y2,'log',stats2,'constant','off');
      Xtemp = repmat(spline2,1,100);
      XX.AAC = Xtemp';

      XX.null = XXC; XX.PAC = XX1;XX.PACAAC = XX3;

      XX.ampAXIS = ampAXIS; XX.phi0 = phi0;
      XX.rpac = max(max((abs(1-XX1./XXC))));
      XX.raac = max(max((abs(1-Xtemp'./XXC))));
      XX.rpacaac = max(max((abs(1-XX3./XXC))));
      [m,I] = max(abs(1-XX.PACAAC./XX.null));
      [M,j] = max(m); %j ampLO, I(j) phiLO
      XX.relPAC = max(abs(XX.PAC(I(j),j)-XX.null(I(j),j)))/max(abs(XX.PACAAC(I(j),j)-XX.null(I(j),j)));
      XX.relAAC = max(abs(XX.AAC(I(j),j)-XX.null(I(j),j)))/max(abs(XX.PACAAC(I(j),j)-XX.null(I(j),j)));

  %Define dense phase points for interpolation.
  phi0 = linspace(-pi,pi,100);
  X0 = spline_phase0(phi0',nCtlPts);
  Amax = max(ampLO); Amin = min(ampLO); stepsize = (Amax-Amin)/99;
  X2eval = Amin:stepsize:Amax; %evaluate on all amplitudes
  X2eval = [ones(size(phi0))',X2eval'];
  X3eval = [X0,AFIX*ones(size(phi0))',AFIX*sin(phi0'),AFIX*cos(phi0')];

  %Determine spline fit and CI.
  [spline0, dylo, dyhi] = glmval(b1,X0,'log',stats1,'constant', 'off');
  splineU = spline0+dyhi;
  splineL = spline0-dylo;
  
  [spline2, dylo2, dyhi2] = glmval(b2,X2eval,'log',stats2,'constant', 'off');
  splineU2 = spline2+dyhi2;
  splineL2 = spline2-dylo2;
  
  [spline3, dylo3, dyhi3] = glmval(b3,X3eval,'log',stats3,'constant', 'off');
  splineU3 = spline3+dyhi3;
  splineL3 = spline3-dylo3;

  %Determine null fit and CI.
  [null0, dylo, dyhi] = glmval(bC,ones(size(phi0)),'log',statsC,'constant', 'off');
  nullU = null0+dyhi;
  nullL = null0-dylo;

  %Find the max absolute percentage change between the two models.
  [r,imx] = max(abs(1-spline0./null0));
  [r2,imx2] = max(abs(1-spline2./null0));
  [r3,imx3] = max(abs(1-spline3./null0));

  %Determine CI for the measure r.
  M = 10000;
  bMC = b1*ones(1,M) + sqrtm(stats1.covb)*normrnd(0,1,nCtlPts,M);
  splineMC = glmval(bMC,X0,'log',stats1,'constant', 'off');
%   nullMC   = mean(splineMC,1);
  mx = zeros(M,1);
  for k=1:M
      %mx(k) = max(abs(1-splineMC(:,k)./nullMC(k)));
      mx(k) = max(abs(1-splineMC(:,k)./splineC));
  end
  r_CI = quantile(mx, [0.025, 0.975]);
  
  %and for r2
  M = 10000;
  bMC = b2*ones(1,M) + sqrtm(stats2.covb)*normrnd(0,1,2,M);
  splineMC = glmval(bMC,X2eval,'log',stats2,'constant', 'off');
  %nullMC   = mean(splineMC,1);
  mx = zeros(M,1);
  for k=1:M
      %mx(k) = max(abs(1-splineMC(:,k)./nullMC(k)));
      mx(k) = max(abs(1-splineMC(:,k)./null0));
  end
  r2_CI = quantile(mx, [0.025, 0.975]);
  
  %and for r3
  M = 10000;
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

  
    if isempty(varargin) || ~strcmp(varargin{1}, 'noplot')
      
      plot(phi0, null0, 'k', 'LineWidth', 2);  hold on
      plot(phi0, nullL, 'k:', 'LineWidth', 2)
      plot(phi0, nullU, 'k:', 'LineWidth', 2)
      
      plot(phi0,spline3,'b','LineWidth',2)
      plot(phi0,splineU3,'b:','LineWidth',2)
      plot(phi0,splineL3,'b:','LineWidth',2)
      
      plot([phi0(imx3) phi0(imx3)], [null0(imx3), spline2(imx3)], 'LineWidth', 2)
      %plot([phi0(imx) phi0(imx)], [null0(imx), spline0(imx)], 'LineWidth', 2)
      hold off
      ylabel('Amplitude')
      xlabel('Phase')
      axis tight
    end

R.PAC = r;
R.AAC = r2;
R.PACAAC = r3;
R.PACci = r_CI;
R.AACci = r2_CI;
R.PACAACci = r3_CI;
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

