function [x1new] = make_pink_noise(alpha,L,dt)

  %alpha=0.33;

  x1 = randn(L,1);
  xf1 = fft(x1);
  A = abs(xf1);
  phase = angle(xf1);

  df = 1.0 / (dt*length(x1));
  faxis = (0:length(x1)/2)*df;
  faxis = [faxis, faxis(end-1:-1:2)];  %(end-1:-1:2)
  oneOverf = 1.0 ./ faxis.^alpha;
  oneOverf(1)=0.0;

  Anew = A.*oneOverf';
  xf1new = Anew .* exp(i*phase);
  x1new = real(ifft(xf1new))';
  
end