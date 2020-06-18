
function cosmic_variance, Cl
vs= size(Cl)
Nl = vs[1]
l = lindgen(Nl) + 2
ll = l*(l+1)
 ll2 = 2.*l+1
sigcv = cl * sqrt(2./ll2)
return, sigcv
end