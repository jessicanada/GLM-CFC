function [mi,p,MI] = modulation_index(Vlo, Vhi,pvals)

  phi=angle(hilbert(Vlo));%Compute phase of low freq signal.
  amp=abs(hilbert(Vhi));	%Compute amplitude of high freq signal.
  
  N      = 18;                        %Number of phase bins.
  p_bins = linspace(-pi,pi,N+1);       %Define the phase bins.
  a_mean = zeros(length(p_bins)-1,1);	%Vector for average amps.
  p_mean = zeros(length(p_bins)-1,1);	%Vector for phase bins.
  
  for k=1:N                           %For each phase bin,
      pL = p_bins(k);					%... lower phase limit,
      pR = p_bins(k+1);				%... upper phase limit.
      indices=find(phi>=pL & phi<pR);	%Find phases falling in bin,
      a_mean(k) = mean(amp(indices));	%... compute mean amplitude,
      p_mean(k) = mean([pL, pR]);		%... save center phase.
  end
  
  %Difference between max and min modulation.
  P = a_mean / sum(a_mean);
  H = -sum(P.*log(P));
  mi= (log(N)-H)/log(N);
  
  p =[];
  if exist('pvals','var') && strcmp(pvals, 'pvals')
      MI = zeros(1,1000);
       for i = 1:1000
             Vhi_prime = AAFT(Vhi,1);
             MI(i) = modulation_index(Vlo,Vhi_prime);%compute MI between Vhi and shifted Vlo
       end                                 %make sure only looking when Vlo shifted
       p = max(.5,length(find(MI>mi)))/length(MI);
  end
  
end