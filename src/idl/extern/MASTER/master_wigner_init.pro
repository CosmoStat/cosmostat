pro master_wigner_init,lmax,lnwpro
k=indgen(2+lmax*3)
lnwpro=lngamma(2*k+1d)-2*lngamma(k+1d)
end
