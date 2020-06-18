pro master_anafast_h,map,ell,cell,lmax,order=order,tmpdir=tmpdir,alms=alms,diranafast=diranafast

if not keyword_set(tmpdir) then tmpdir='/tmp/'

truc=strtrim(long(randomu(seed,1)*1e6),2)
fitsmap=tmpdir+'/tmpmap'+truc(0)+'.fits'
fitscl=tmpdir+'/tmpcl'+truc(0)+'.fits'
txtfile=tmpdir+'/txt'+truc(0)+'.txt'
if (keyword_set(alms)) then almsfile = tmpdir+'/alms'+truc(0)+'.txt'

spawn,'rm -f '+fitsmap+' '+fitscl+' '+txtfile

if not keyword_set(order) then order=1
if order eq -1 then order=0

nside=npix2nside(n_elements(map))

;if lmax gt 3*nside-1 then lmax=3*nside-1

;add_nside_fits,info_hdr,nside=nside,error=error
;nside_header=['NSIDE   =                  '+strtrim(nside,2)+' /Binary table']
write_fits_map,fitsmap,map,/nested
;temp=mrdfits(fitsmap,1,h)
;sxaddpar,h,"NSIDE",nside
;write_fits_map,fitsmap,map,h,/ring
order = 0
;close,1
;openw,1,txtfile
;printf,1,'1'
;printf,1,fitsmap
;printf,1,lmax
;printf,1,"''"; file mask1
;printf,1,'0.'
;printf,1,order
;printf,1,"''" ;precomputed plms file
;printf,1,fitscl
;if (not keyword_set(alms)) then printf,1,' ' else printf,1,almsfile
;printf,1,'2'
;close,1

;close,1
;openw,1,txtfile
;printf,1,'1'

;printf,1,'infile='+fitsmap
;printf,1,'polarisation=false'
;printf,1,'nlmax='+string(lmax)
;printf,1,"''"; file mask1
;printf,1,'iter_order=2'
;printf,1,'0.'
;printf,1,order
;;printf,1,"''" ;precomputed plms file
;printf,1,'outfile='+fitscl
;if (not keyword_set(alms)) then printf,1,' ' else printf,1,almsfile
;printf,1,'won=2'
;close,1

close,1
openw,1,txtfile
printf,1,'infile='+fitsmap
printf,1,'polarisation=false'
printf,1,'nlmax='+string(lmax)
printf,1,'iter_order=2'
printf,1,'outfile='+fitscl
;if (not keyword_set(alms)) then printf,1,' ' else printf,1,almsfile
printf,1,'won=2'
printf,1,'double_precision=true'

close,1



spawn, diranafast+'anafast_cxx '+txtfile
cl = mrdfits(fitscl,1,temphh)
cell = cl.(1)

spawn,'cp -f  '+fitscl+' cl2_cxx.fits'
spawn,'rm -f '+fitsmap+' '+fitscl+' '+txtfile
if keyword_set(alms) then begin 
  alms=mrdfits(almsfile,1)
  spawn,'rm -f '+almsfile
endif

ell=dindgen(n_elements(cell))

end
