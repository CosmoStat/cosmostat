function master_make_mll,ell,well,lmax

ellm=ell(0:lmax)
master_wigner_init,n_elements(ell),lnwpro
m=dblarr(n_elements(ellm),n_elements(ellm))
;plot,findgen(n_elements(ellm)),findgen(n_elements(ell)),/nodata, $
;  title='Computing the Mll matrix ...',/xsty,/ysty
for l1=0,n_elements(ellm)-1 do begin & $
;  plots,ellm(l1),ellm(l1),psym=1 & $
  c=(2*ell+1)*well & $
  for l2=0,n_elements(ellm)-1 do begin & $
     b=c*(master_wigner3j2(ellm(l1),ellm(l2),ell,lnwpro)) & $
     a=total(b) & $
     m(l1,l2)=(2*ellm(l2)+1)*a/(4.*!pi) & $
  end & $
end
return,m
end
