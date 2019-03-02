function [Vlo, Vhi, t, Alo] = generate_interesting_vlo_and_vhi(CFC_parameter)

  N  = 60000+4000;                        %Number of data poinsts to use.
  dt = 0.002;  Fs = 1/dt;  fNQ=Fs/2;      %Sampling parameters
  
  % ------------------------------------- Make Vlo.
  Vpink = make_pink_noise(1,N,dt);
  Vpink = Vpink - mean(Vpink);
  %Filter into low freq band (4-7 Hz)
  locutoff = 4;
  hicutoff = 7;
  filtorder = 3*fix(Fs/locutoff);
  MINFREQ = 0;
  trans          = 0.15;
  f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
  m=[0       0                      1            1            0                      0];
  filtwts = firls(filtorder,f,m);
  Vlo = filtfilt(filtwts,1,Vpink);
  
  % ------------------------------------- Make Vhi.
  Vpink = make_pink_noise(1,N,dt);
  Vpink = Vpink - mean(Vpink);
  %Filter into high frequency band (100-140 Hz)
  locutoff = 100;
  hicutoff = 140;
  filtorder = 10*fix(Fs/locutoff);
  MINFREQ = 0;
  trans          = 0.15;
  f=[MINFREQ (1-trans)*locutoff/fNQ locutoff/fNQ hicutoff/fNQ (1+trans)*hicutoff/fNQ 1];
  m=[0       0                      1            1            0                      0];
  filtwts = firls(filtorder,f,m);
  Vhi = filtfilt(filtwts,1,Vpink);
  
  % ------------------------------------- Drop the edges to remove filter artifacts.
  Vlo = Vlo(2001:end-2000);
  Vhi = Vhi(2001:end-2000);
  t   = (1:length(Vlo))*dt;
  
  % ------------------------------------- Create a signal with clever PAC.
  
  %Define a signal that is "on" and "off" at specific Vlo phases.
  %s = cos(angle(hilbert(Vlo))+pi/4);
  
  %Define a signal that has only positive (no negative) relationships
  [~, ipks] = findpeaks(Vlo);
  threshold = 0;
    AmpLo = abs(hilbert(Vlo));
    s = zeros(size(Vhi));                               % Define empty modulation envelope.
    for i0=1:length(ipks)                               % At every low freq peak,
        if ipks(i0) > 10 && ipks(i0) < length(Vhi)-10   % ... if indices are in range of vector length.
%             val = AmpLo(ipks(i0));
%             if val>=quantile(AmpLo(ipks),threshold) %If the amplitude is in the top half of the distribution
%                 scalar = 1;               %Make the modulation large
%             else
%                 scalar = (val/max(Vlo))^6; %Otherwise, make the modulation small (cube of small value<1)
%                 %scalar = 0;
%             end
%             s(ipks(i0)-10:ipks(i0)+10) = scalar*hann(21); %Scaled Modulation
              s(ipks(i0)-10:ipks(i0)+10) = 1*hann(21); %Scaled Modulation
        end
    end
    s = s/max(s); %normalize so it falls between 0 and 1
   
  
  %Get the Vlo amplitude.
  Alo = abs(hilbert(Vlo));
  
  %Create a matrix of predictors.
  %      Predictor #1: constant
  %      Predictor #2: Alo * s = this is an AmpLo * PhaseLo term.
  
  X   = [ones(length(Alo),1), (Alo'-mean(Alo)).*s'];
  %NOTE:                      the 2nd term is Alo minus its mean!
  
  %Create Ahi.
  impact_Alo_and_PhiLo_on_Ahi = CFC_parameter;	% Adjust this to create PAC. Try = 0 vs = 5.
  b = [1, impact_Alo_and_PhiLo_on_Ahi]';                  % Define coefficeints of predictor
  Ahi = glmval(b,X,'log','constant','off');               % Generate Ahi via a GLM.
  Vhi = 0.001* Ahi .* cos(angle(hilbert(Vhi')));          % ...  create Vhi, using Ahi we just created and original Vhi phase.

end