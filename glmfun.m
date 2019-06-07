function [XX,P] = glmfun(Vlo, Vhi,pval,ci,AIC,varargin)
%INPUTS:
% Vlo:                   Low frequency signal
% Vhi:                   High frequency signal
% pval:                  'pval' gives bootstrapped p-values for R
%                        'none' outputs no p-values (faster)
% ci:                    'ci' gives confidence intervals, 'none' gives no confidence intervals
% AIC:                   'AIC' performs an AIC-based selection procedure to
%                        choose the number of splines
%                        'none' chooses the standard (n=10) number of splines
%                        (faster)
% varargin:              optionally, include the parameter q indicating which quantiles
%                        of AmpLo you'd like to fit over

%OUTPUTS:
% XX.rpac:          R_PAC value, confidence intervals XX.rPAC_CI
% XX.raac:          R_AAC value, confidence intervals XX.rAAC_CI
% XX.philow:        3D surface for Phi_low model in Phi_low, A_low, A_high space
% XX.alow:          3D surface for A_low model in Phi_low, A_low, A_high space
% XX.alowphilow:    3D surface for A_low,Phi_low model in Phi_low, A_low, A_high space
% P.rpac:           p-value for RPAC statistic
% P.raac:           p-value for RAAC statistic
      
  nCtlPts = 10;

  %Compute phase and amplitude.
  phi = angle(hilbert(Vlo));
  amp = abs(hilbert(Vhi));
  ampLO = abs(hilbert(Vlo));
  
  if exist('AIC','var') && strcmp(AIC, 'AIC')

    Y = amp';
    CtlPts = (4:1:30);
    %Compute the AIC.
    AIC1 = zeros(size(CtlPts)); %AIC3 = zeros(size(CtlPts));
    for k=1:length(CtlPts) %for each suggested # knots
        nCtlPts = CtlPts(k);
        X1 = spline_phase0(phi',nCtlPts);

        [~,  dev1, ~]  = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off'); %PAC

        AIC1(k) = dev1 + 2*nCtlPts; %compute AIC
    end
    %Select the # control points from AIC, and plot the AIC.
    [~, imn1] = min(AIC1);

    nCtlPts1 = CtlPts(imn1);

    figure(1); clf();
    plot(CtlPts, AIC1, 'k', 'LineWidth',2)
    hold on
    plot([nCtlPts1, nCtlPts1], [min(AIC1) max(AIC1)], 'r', 'LineWidth', 2)
    hold off
    axis tight
    xlabel('# control points')
    ylabel('AIC')
    set(gca,'FontSize',13)
% 
     fprintf(['Suggested number of control points  is ' num2str(nCtlPts1) '\n'])

    nCtlPts = nCtlPts1;
  end
  
  %Define variables for GLM procedure.
  Y = amp';                                                                 %high frequency amplitude
  X1 = spline_phase0(phi',nCtlPts);                                         %low frequency phase
  X2 = [ones(size(Y)),ampLO'];                                              %low frequency amplitude
  X3 = [X1,ampLO',sin(phi').*ampLO',cos(phi').*ampLO'];                     %low frequency phase, low frequency amplitude, interaction terms
  XC = ones(size(Y));                                                       %ones (null)

  %Perform GLM.
  [b1,  dev1, stats1]  = glmfit(X1,  Y, 'gamma', 'link', 'log', 'constant', 'off'); %PAC
  [b2,  dev2, stats2] = glmfit(X2, Y, 'gamma','link','log','constant','off');       %AAC
  [b3,  dev3, stats3] = glmfit(X3, Y, 'gamma','link','log','constant','off');       %CFC
  [bC, dev0, statsC] = glmfit(XC, Y, 'gamma', 'link', 'log', 'constant', 'off');    %null
  
  
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
      Y2 = [ones(size(temp')),temp'];
      [spline2,~,~] = glmval(b2,Y2,'log',stats2,'constant','off'); %AAC model, function of Alow
      Xtemp = repmat(spline2,1,100);                                         %AAC model constant in PhiLow dimension
      
      XX.alow = Xtemp'; XX.null = XXC; XX.philow = XX1;XX.alowphilow = XX3;          %3D model surfaces
         
      XX.ampAXIS = ampAXIS; XX.phi0 = phi0;     %axes

      XX.rpac = max(max(abs(1-Xtemp'./XX3)));

      XX.raac = max(max((abs(1-XX1./XX3))));

  
  if exist('pval','var') && strcmp(pval, 'pval')
    M = minvals(Vlo,Vhi); %find empirical p-values
    P.rpac = max(.5,length(find(M.rpac>XX.rpac)))/length(M.rpac);
    P.raac = max(.5,length(find(M.raac>XX.raac)))/length(M.raac);
  else
    P = 'No p-values output';
  end
  
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
    spline_a = glmval(bMC,Y_a,'log',stats3,'constant','off');
    
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
function M = minvals(Vlo,Vhi)
    K = 100;
    RPAC = zeros(1,K); RAAC = zeros(1,K);
    N = zeros(1,K); L = zeros(1,K);
    for i = 1:K
        Vhi_prime = AAFT(Vhi,1);
        [XX] = glmfun(Vlo,Vhi_prime','none','none','none',.05);        %compute R statistics between Vhi and shifted Vlo
        RPAC(i) = XX.rpac;
        RAAC(i) = XX.raac;
    end
    M.rpac = RPAC;  M.raac = RAAC; M.shiftN = N; M.shiftL = L;
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

