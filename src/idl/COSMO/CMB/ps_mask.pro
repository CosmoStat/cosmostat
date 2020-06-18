pro ps_mask,freq,Slim,nside=nside,gb=gb,n_sigma_T=n_sigma_T,n_sigma_P=n_sigma_P,fwhm=fwhm
;+
; NAME:
;    ps_mask
; PURPOSE:
;    Computes a mask map of PS brigther than Slim at a given frequency
; CALLING SEQUENCE:
;    ps_mask, freq, Slim,[nside=nside, gb=gb, n_sigma=n_sigma, n_sigma_P=n_sigma_P, fwhm=fwhm]
; INPUT PARAMETERS:
;    freq  = frequency at which calculation is to be done [GHz]
;    Slim  = PS brighter than this flux limits are masked [mJy]
; OPTIONAL INPUT PARAMETERS
;	 nside=  To use a different nside of the choosen frequency
;        gb = To consider only PS at /Glat/>gb [degree]
;        n_sigma_T, n_sigma_P = If set, Slim is ignored and recomputed
;        corresponding to n_sigma_T * sigma_T or n_sigma_P * sigma_P
;        where sigma_T or sigma_P is the cmb rms at frequency=freq in
;        total or polarized intensity, respectively;
;        the cmb sigma_T and sigma_P at nside=2048 are those of the Planck
;        reference sky emission available at www.sissa.it/~planck/reference_sky/CMB/;
;        at nside<2048, they are evaluated using the HEALPix ud_grade routine
;	 fwhm= To use a different fwhm of the choosen frequency [arcmin]
; OUTPUT:
;       A FITS file with a Mask (pixels around 2 beam-size of a PS
;       brigter than Slim has value 0 and 1 otherwise) at the choosen frequency
; EXAMPLES:
;     Mask at 70 GHz with PS>200 mJy
;
;     	        ps_mask,70.,200.
;
;     Mask at 217 GHz with PS>500 mJy and /Glat/>5deg
;
;     	        ps_mask,217.,500.,gb=5
;
;     Mask at 30 GHz for S>1000 mJy and nside 2048
;
;     	        ps_mask,30.,1000.,nside=2048.
;
;     Mask at 100 GHz for cut at 5 times the sigma of CMB and nside 2048
;
;     	        ps_mask,30.,1000.,nside=2048.,n_sigma=5
;
;     Mask at 217 GHz for S>1000 mJy and fwhm=7.1'
;
;     	        ps_mask,217.,1000.,fwhm=7.1
;
; REVISION HISTORY:
;    July 28, 2006, edited to match the format of the working group 2
;    CMB and foreground template repository, prepared by Joaquin
;    Gonzalez-Nuevo
;    Agoust 8, 2006, n_sigma optional input added by Carlo Baccigalupi
;    Septembre 1, 2006, fwhm optional input added by Joaquin Gonzalez-Nuevo
;    September 6, 2006, fwhm info added in the output file name by F. Stivoli
;    October 3, 2006, BUG corrected: ARs is not in seconds but in deci-seconds
;-

;freq_index
nu=[30,44,70,100,143,217,353,545,857]
kk=min(abs(nu-freq),id_nu)

; Beam
if defined(fwhm) then begin
	beam=fwhm/60./180.*!PI/sqrt(2.*alog(2.))/2.
endif else begin
	fwhm_nu=[33.,24.,14.,10.,7.1,5.,5.,5.,5.]
	beam=fwhm_nu[id_nu]/60./180.*!PI/sqrt(2.*alog(2.))/2.
        fwhm=fix(fwhm_nu[id_nu])    ;used later to create the file name
endelse
print, beam

;nside
if defined(nside) then begin
	nside=nside
endif else begin
	nside_nu=[256,512,1024,1024,2048,2048,2048,2048,2048]
	nside=nside_nu[id_nu]
endelse

;n_sigma
; units and useful quantities for the conversion from micro K thermodynamic to mJy/sr:
h=6.626176e-27                  ; erg * s
k=1.380662e-16                  ; erg / K
c=2.99792458e10                 ; cm / s
T_cmb=2.726                     ; K
conversion_microKthermodynamic_K=1.e-6*(h*freq*1.e9/k/T_cmb)^2*exp(h*freq*1.e9/k/T_cmb)/(exp(h*freq*1.e9/k/T_cmb)-1.)^2
conversion_K_mJy_sr=2.*k*(freq*1.e9/c)^2*1.e26
conversion_microKthermodynamic_mJy_sr=conversion_microKthermodynamic_K*conversion_K_mJy_sr
if (nside eq 2048) then begin
sigma_P=6.29826			; micro-K thermodynamic temperature
sigma_T=110.176			; micro-K thermodynamic temperature
endif
if (nside eq 1024) then	begin
sigma_P=6.08849			; micro-K thermodynamic temperature
sigma_T=109.552			; micro-K thermodynamic temperature
endif
if (nside eq 512) then begin
sigma_P=5.43045			; micro-K thermodynamic temperature
sigma_T=107.475			; micro-K thermodynamic temperature
endif
if (nside eq 256) then begin
sigma_P=4.06417			; micro-K thermodynamic temperature
sigma_T=102.282			; micro-K thermodynamic temperature
endif
if (nside eq 128) then begin
sigma_P=2.42632			; micro-K thermodynamic temperature
sigma_T=92.9725			; micro-K thermodynamic temperature
endif
if (nside eq 64) then begin
sigma_P=1.14274			; micro-K thermodynamic temperature
sigma_T=78.5644			; micro-K thermodynamic temperature
endif
if (nside eq 32) then begin
sigma_P=0.679663		; micro-K thermodynamic temperature
sigma_T=62.3515			; micro-K thermodynamic temperature
endif
if (nside eq 16) then begin
sigma_P=0.435442		; micro-K thermodynamic temperature
sigma_T=51.1995			; micro-K thermodynamic temperature
endif
if (nside eq 8) then begin
sigma_P=0.365490		; micro-K thermodynamic temperature
sigma_T=42.5380			; micro-K thermodynamic temperature
endif
if (nside eq 4) then begin
sigma_P=0.326714		; micro-K thermodynamic temperature
sigma_T=33.7957			; micro-K thermodynamic temperature
endif 
if defined(n_sigma_T) then begin
        Slim=n_sigma_T*sigma_T*conversion_microKthermodynamic_mJy_sr*4.*!pi/12./nside/nside
print,'flux limit corresponding to',n_sigma_T,'sigma for total intensity cmb at',freq,'GHz:',Slim,'mJy'
endif
if defined(n_sigma_P) then begin
        Slim=n_sigma_P*sigma_P*conversion_microKthermodynamic_mJy_sr*4.*!pi/12./nside/nside
print,'flux limit corresponding to',n_sigma_P,'sigma for polarized intensity cmb at',freq,'GHz:',Slim,'mJy'
endif 

;reading catalogues
cat_radio=readfits('radiosource_catalog.fits')
;cat_radio[0,*]=cat_radio[0,*]/24.*360. ;hms-> deg
cat_radio[2,*]=cat_radio[2,*]/1000. ;mJy-> Jy

cat_IR=fltarr(500000L,3)

i=0L
sig='-'
close,1
openr,1,'IRAS_PSC.DAT'
while not eof(1) do begin
	readf,1,arh,arm,ards,sig,decd,decm,decs,S12,S25,S60,S100,format='(11x,I2,I2,I3,A1,I2,I2,I2,11x,E9.3,E9.3,E9.3,E9.3)'
	sign=+1
	if(sig eq '-') then sign=-1.
	cat_IR[i,*]=[arh+arm/60.+ards/3600./10.,sign*(decd+decm/60.+decs/3600.),S60]
	i=i+1
endwhile
close,1
cat_IR=cat_IR[0:i-1,*]

;(AR,Dec) to  (GLon,GLat)
call_procedure,'glactc',cat_radio[0,*],cat_radio[1,*],2000,glon,glat,1
cat_radio[0,*]=glon
cat_radio[1,*]=glat
glon=0. & glat=0.
call_procedure,'glactc',cat_IR[*,0],cat_IR[*,1],1950,glon,glat,1
cat_IR[*,0]=glon
cat_IR[*,1]=glat
glon=0. & glat=0.


;<S/S60>=[100GHz,143,217,353,545,857]
S_S60=[0.0002,0.0008,0.0041,0.0241,0.1071,0.4229]
cat_IR[*,2]=cat_IR[*,2]*S_S60[(id_nu-3)>0]

;Flux selection
w_r=where(cat_radio[2,*] ge Slim/1000.)
cat_radio=cat_radio[*,w_r]
w_ir=where(cat_IR[*,2] ge Slim/1000.)
cat_IR=cat_IR[w_ir,*]
w_r=0. & w_ir=0.

;GLat selection
if defined(gb) then begin
	w_r=where(abs(cat_radio[1,*]) ge gb)
	cat_radio=cat_radio[*,w_r]
	w_ir=where(abs(cat_IR[*,1]) ge gb)
	cat_IR=cat_IR[w_ir,*]
	w_r=0. & w_ir=0.
endif

;Healpix pixels
ang2pix_ring, nside, (90.-cat_radio[1,*])/180.*!PI, cat_radio[0,*]/180.*!PI, ipring_r
ang2pix_ring, nside, (90.-cat_IR[*,1])/180.*!PI, cat_IR[*,0]/180.*!PI, ipring_ir

Map=fltarr(12.*nside*nside)

Map[*]=1.

;Healpix circles
for i=0L,n_elements(ipring_r)-1 do begin
	pix2vec_ring,nside,ipring_r[i],vector
	query_disc,nside,vector,2*beam,listpix_r
	Map[listpix_r]=0.
endfor
for i=0L,n_elements(ipring_ir)-1 do begin
	pix2vec_ring,nside,ipring_ir[i],vector
	query_disc,nside,vector,2*beam,listpix_ir
	Map[listpix_ir]=0.
endfor

; Sigma estimation
counts=dblarr(100,10)

i=0L
close,1
openr,1,'c_dezotti_all.dat'
while not eof(1) do begin
	readf,1,h,v1,v2,v3,v4,v5,v6,v7,v8,v9,format='(F5.1,1x,F7.4,1x,F7.4,1x,F7.4,1x,F7.4,1x,F7.4,1x,F7.4,1x,F7.4,1x,F7.4,1x,F7.4)'
	counts[i,*]=[h,v1,v2,v3,v4,v5,v6,v7,v8,v9]
	i=i+1
endwhile
close,1
counts=counts[0:i-1,*]
w=where(counts[*,0] gt -6 and counts[*,0] lt alog10(Slim/1000.))
I=total(10^counts[w,id_nu+1]*(10^counts[w,0])^2)
sigma_radio=sqrt(I*(4*!PI/(12.*nside*nside)))

print,'*************************************************************'
print,'Point Source Mask at',freq,' GHz'
print,'nside:',nside,', Slim:',Slim,' mJy'
if defined(gb) then print,'galaxy cut:',gb,' deg'
print,'Radio Sources > Slim:',n_elements(ipring_r)
print,'IR Sources > Slim:',n_elements(ipring_ir)
npix=n_elements(where(Map lt 1.))
print,'Masked pixels:',npix,'(',float(npix)/float(n_elements(Map))*100.,'%)'
if defined(gb) then print,'Masked pixels fraction outside the cut:','(',float(npix)/float(n_elements(Map))/(1.-cos(gb*!pi/180.))*100.,'%)'
print,'Sigma of undetected RADIO sources:',sigma_radio*1000.,' mJy'
print,'*************************************************************'

outputfile='./mask_ps_'+strtrim(round(freq),2)+'GHz_Slim_'+strtrim(round(Slim),2)+'mJy_beam_'+strtrim(fwhm,2)+'amin_nside'+strtrim(nside,2)+'.fits'
print,'Writing to ',outputfile

vs = size(Map)
nx = vs[1]
ir = lindgen(nx)
RING2NEST, Nside, ir, in
mollview, map, /on
Mapn = Map 
Mapn[in] = Map
mrs_tv, mapn
write_fits_map,outputfile,Map,/ring,coord='G'

;mollview,Map,/online
end
