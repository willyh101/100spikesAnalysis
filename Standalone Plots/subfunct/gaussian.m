function val = gaussian(x,amp,sig)
fudge = 1e-3;
s2 = sig^2 + fudge;
val = amp*exp(-0.5*x.^2/s2)./sqrt(2*pi*s2);
end