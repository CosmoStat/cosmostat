;+
; NAME: 
;       IM_ISOSPEC
;
; PURPOSE: 
;       Calculate the isotropic power spectrum.
;       
; CALLING:
;
;       IM_ISOSPEC, Data, powspec, x=x, l=l, fit=fit, tit=tit
;
; INPUTS:
;       Data: image to analyze
;
; INPUT KEYWORDS:
;           fit: float or fltarr(2): the log(power_spectrum) is fitted 
;                                    by a line between [fit:*] or between
;                                    [fit(0):fit(1)]
;
;           plot: the log power spectrum is plotted
;
;           l: integer: if set, the power spectrum is multiplied by l(l+1)
;                        where l is the frequency.
;
; OUTPUT KEYWORDS:
;           x: 1D IDL array: x-axis (frequencies)
;        
; OUTPUTS:
;           powspec: 1D IDL array: power spectrum
;
; EXTERNAL CALLS:
;           im_isospec (C++ program)
;
; EXAMPLE:
;      im_isospec, Image, powspec, /plot, fit=0.2
;      im_isospec, Image, powspec, /plot, fit=[0.2,0.5]
;
;-



;===================================================================================================
;
; NAME: 
; 				POLAR_FFT
;
;
; PURPOSE: 
; 				process the polar fft of an image
;
; CALLING: 
; 				polar_fft, f, m, real, ima
;
; INPUTS: 
;					f --- input image
; 				m --- cut-off parameter in the approximation
; 				(the larger it is the better approximation is)
;	
; OUTPUT: 
; 				real --- real part of the polar fft
; 				ima --- imaginary part of the polar fft 
;					(each line of the polar fft is a radius, 
; 				the line 0 corresponds to angle -!pi/2)
;
; HISTORY:
;	Written: Sandrine Pires Jan 2007.
;-
;-------------------------------------------------------------------------------

function gettmpfilename
; spawn,'mktemp',filename_tmp
if  !version.os EQ 'linux' then spawn,'mktemp -p /dev/shm', filename_tmp $
else  spawn,'mktemp',filename_tmp 

; filename_tmp = strcompress('tmp_in_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
return, filename_tmp
end

;===============================================================


pro delete, filename
case !version.os of
    'vms': cmd = 'delete'
    'windows': cmd = 'del'
    'MacOS': goto, notsup
    else: cmd = '\rm -f'
    endcase
cmd = cmd + ' ' + filename  
spawn, cmd  
goto, DONE
NOTSUP: print, "This operation is not supported on the Macintosh'
DONE:
end




;===================================================================================================
;
; NAME: DAT2FITS
;
;
; PURPOSE: to change a .dat file (from matlab save) in a .fits file
;
;
; CALLING: dat2fits, file_dat=file_dat filefits file_fits=file_file
;
; INPUTS: 
;	
;	
; OUTPUT: 
;	
;
; HISTORY:
;	Written: Sandrine Pires Jan 2007.
;-
;-------------------------------------------------------------------------------

pro dat2fits, file_dat=file_dat, filefits, file_fits=file_fits

tab = read_ascii(file_dat)

filefits = tab.(0)
 
writefits, file_fits, filefits

end

;===================================================================================================
; NAME:FITS2DAT
;
;
; PURPOSE: to change a .dat file (from matlab save) in a .fits file
;
;
; CALLING: 
; 					fits2dat, file_fits=file_fits, filefits, file_dat=file_dat
; INPUTS: 
;	
;	
; OUTPUT: 
;	
;
; HISTORY:
;	Written: Sandrine Pires Jan 2007.
;-
;-------------------------------------------------------------------------------

pro fits2dat, file_fits=file_fits, filefits, file_dat=file_dat

tab = readfits(file_fits)
openw, 1, file_dat
printf, 1, double(tab)
close, 1
free_lun, 1

end

;===============================================================

pro polar_fft, f, real, ima, PolarFFT_CutOffParam=PolarFFT_CutOffParam

if keyword_set(PolarFFT_CutOffParam) then  m = PolarFFT_CutOffParam else m = 8

;a = systime(1)
chemin = './'
sz =size(f)
N = sz[1]
T=3*N; %T=5*N/2
R=3*N/2

fr=float(f) 
fi=imaginary(f)

IN_FITSFile_RE = gettmpfilename()
IN_FITSFile_IM = gettmpfilename()
IN_DATFILE_RE = 'input_data_r.dat'
IN_DATFILE_IM = 'input_data_i.dat'
OUT_DAT_FILE_RE = 'polar_fft_r.dat'
OUT_DAT_FILE_IM = 'polar_fft_i.dat'
OUT_FITS_FILE_RE = gettmpfilename()
OUT_FITS_FILE_IM = gettmpfilename()

writefits,  IN_FITSFile_RE, fr
fits2dat, file_fits=IN_FITSFile_RE , filefits, file_dat=IN_DATFILE_RE 
writefits, IN_FITSFile_IM, fi
fits2dat, file_fits=IN_FITSFile_IM, filefits, file_dat=IN_DATFILE_IM

com = 'cea_polar_fft' +' ' +string(N)+' ' +string(T)+' ' +string(R)+' ' +string(m)
spawn, com
;dat2fits, file_dat='polar_fft_error.dat',filefits,file_fits='polar_fft_error.fits'
;polar_fft_error = readfits('polar_fft_error.fits')

dat2fits, file_dat=OUT_DAT_FILE_RE, filefits, file_fits= OUT_FITS_FILE_RE 
dat2fits, file_dat= OUT_DAT_FILE_IM ,filefits,file_fits= OUT_FITS_FILE_IM

real = readfits(OUT_FITS_FILE_RE)
ima = readfits(OUT_FITS_FILE_IM )
delete,  IN_FITSFile_RE
delete,  IN_FITSFile_IM
delete,  IN_DATFILE_RE
delete,  IN_DATFILE_IM
delete,  OUT_DAT_FILE_RE
delete,  OUT_DAT_FILE_IM
delete,  OUT_FITS_FILE_RE
delete,  OUT_FITS_FILE_IM


;print, 'time =', systime(1)-a
end


;==========================================================================================

function convert_ps2cl, psize_arcmin, N 
;  CL = convert_ps2cl(psize, N)  * PS

; Size of the field in degree
field_deg = float(N) * psize_arcmin / 60.
field_rad = (field_deg*!pi)/180.

; Surface of the field in radian^2
surf_field_rad = field_rad^2.

;cl normalisation that verify : \sigma = \integral{l*cl/2*!pi} dl
Norm = (surf_field_rad)/float(N)^2
return, Norm
end

;==========================================================================================

function convert_l2x, psize_arcmin, N
; l = convert_l2x(psize, N)  * x 

; Size of the field in degree
field_deg = float(N) * psize_arcmin / 60.
field_rad = (field_deg*!pi)/180.

NormL = (2.*!pi*float(N))/ float(field_rad)

return, NormL
end

;==========================================================================================

function ps2cl, x, powspec, indl=indl, psize_arcmin=psize_arcmin, N=N

if not keyword_set(N) then N = long( max(x+0.5))
if not keyword_set(psize_arcmin) then psize_arcmin =  1.71774   ; pixel_size(2048)

indl = x* convert_l2x(psize_arcmin, N)
;cl normalisation that verify : \sigma = \integral{l*cl/2*!pi} dl
cl = powspec * convert_ps2cl(psize_arcmin, N)

return, Cl
end
;==========================================================================================

function cl2ps, indl, Cl, x=x, psize_arcmin=psize_arcmin, N=N

if not keyword_set(N) then N = long( max(l+0.5))
if not keyword_set(psize_arcmin) then psize_arcmin =  1.71774   ; pixel_size(2048)

x = indl / convert_l2x(psize_arcmin, N)
powspec = Cl / convert_ps2cl(psize_arcmin, N)

return, powspec
end

;==========================================================================================

function hanning, N, WPar=Wpar
  if not keyword_set(WPar) then WPAR = 0.5
  if N mod 2 EQ 1 then N2 = float(N)/2. $
  else N2 = float(N+1)/2.
  
  C = !PI/Wpar;
  WinTab = fltarr(N)
   for i=0,N-1 do begin
      u = (float(i) - N2) / N
      if (ABS(u) < Wpar) then WinTab[i] = 0.5 + 0.5 * cos(C*u) $
      else WinTab[i] = 0;
   end
   return, WinTab
end

;=====================================================================

function win_spline, N
Win = dblarr(N,N)
for i=0L,N-1 do begin
for j=0L,N-1 do begin
  Res = sqrt((i-N/2)^2 + (j-N/2)^2)
  res = 2.0d * res / (N/2)
  Win[i,j] = 1.0d /12.0 * (( abs(res-2d))^3. - 4.0* (abs(res-1))^3 + 6 *(abs(res))^3 - 4.0 *( abs(res+1))^3+(abs(res+2))^3)
  end
end
  ind = where(Win LT 0, c)
  if c GT 0 then Win[ind] = 0
  return, Win / max(Win)
end

;=====================================================================
 
PRO im_isospec, imag, powspec, x=x, plot=plot, l=l, fit=fit,a=a,b=b, tit=tit, nolog=nolog, opt=opt, makefit=makefit, polarfft=polarfft, psize_arcmin=psize_arcmin, Cl=Cl, indl=indl, NormL=NormL, ind2l= ind2l, winAppo=winAppo


if N_PARAMS() LT 2 then begin 
        spawn, 'im_isospec'
        print, 'CALL SEQUENCE: im_isospec, imag, powspec, x=x, plot=plot, fit=fit, l=l'
        goto, DONE
        end

if not keyword_set(Opt) then Opt = ' -r0.5 '  

Nl = (size(imag))[2]
Nc = (size(imag))[1]
 
if Nl eq Nc then begin
  N = nc
endif else begin
  print, "ERROR: Nx non equal to ny"
  goto, DONE
endelse

if keyword_set(winAppo) then begin
   Win = win_spline(N) ; hanning(N)
   Win = Win#Win
   NormValWin = total(win^2) / float(N)^2
end else begin
  NormValWin = 1.
  Win = 1.
end

; By default the pixel size is the PLANCK pixel size (nside=2048) in arc minute
if not keyword_set(psize_arcmin) then psize_arcmin =  1.71774   ; pixel_size(2048)

 
if not keyword_set(polarfft) then BEGIN
   NameImag = 'xx_imag.fits'
   NameResult = 'xx_result.fits'
   writefits, NameImag, imag*Win
   com = 'im_isospec ' + OPT + ' ' + NameImag + ' ' + NameResult
   spawn, com
   Result = readfits(NameResult, /silent)
   x = Result(*,0)
   np = (size(x))[1]
   Result = Result(0:np-2,*)
   x = Result(*,0)
   np = (size(x))[1]

   if keyword_set(l) then powspec = Result(*,2)   $
   else powspec = Result(*,1)  
   delete, NameImag
   delete, NameResult
end else begin ; POLAR FFT
  ;polar fft transform
   polar_fft, imag*Win, real, ima, PolarFFT_CutOffParam=PolarFFT_CutOffParam
   sz =size(real)
   R = sz[1]
   T = sz[2]
   N = T/3 
   spec = fltarr(R/2 + 1)
   som = fltarr(R/2 + 1)
   for i = 0, R-2 do begin
    for j = 0, T-2 do begin
	   ii = abs(i - R/2)
	   som(ii) = som(ii) + 1
     spec(ii) = spec(ii) + real(i,j)^2 + ima(i,j)^2
  endfor
endfor 
powspec = spec/(float(som)*N^2)
x=indgen(R/2 + 1)/ float(R)
end

powspec = powspec / NormValWin

; Compute the Cl from the power spectrum
;l normalisation
indl = x* convert_l2x(psize_arcmin, N)
ind2l = indl*(indl +1)/(2*!pi)
l = indl
cl = powspec * convert_ps2cl(psize_arcmin, N)

if keyword_set(MakeFit) then BEGIN
  F1 = x(1)
  F2 = x[np-1]
  if keyword_set(fit) then begin
     vs = (size(fit))[0]
     if vs EQ 0 then F1 = fit $
     else if vs GE 1 then begin
            F1 = fit[0]
            F2 = fit[1]
          end
   end
     
  ind = where ( x GE F1 and x LE F2, c)
  res = LINFIT( alog(x(ind)), alog(powspec(ind))) 
  a = res(1)
  b =  res(0)
  print, 'LINFIT(Log(PS)=a*Log(freq)+b): a = ', res(1),  ' b = ', res(0)
  
  xx = x
  yy = a*alog(xx) + b

  if keyword_set(plot) then begin
  	 ; plot, alog(x), alog(powspec), xtitle='Log Freq', ytitle='Log PS', xrange=[alog(F1),(F2)]
 	  ; oplot,alog(xx),yy

	  Model = exp(b)*x^a
	  if not keyword_set(tit) then tit=0
 	 if not keyword_set(nolog) then plot, x, alog(powspec), xtitle='Freq', ytitle='Log PS', title=tit $
 	 else plot, x, powspec, xtitle='Freq', ytitle='PS', title=tit

 	  oplot, x, alog(Model), line=2

 	 legend = "a=" + strcompress(string(a,'$(f6.2)'),/remove_all) +  $
          ",b=" + strcompress(string(b,'$(f6.2)'),/remove_all)
 	 if keyword_set(plot) then xyouts,0.5,0.8,legend,/normal

 	  M = max(alog(powspec))
 	  x(*) = 0.25
 	  y=x
 	 y(*) = findgen(np) * (2*M) - M 
 	 ;oplot, x,y
	  x = x / 2
	  ;oplot, x,y
	  x = x / 2
	  ;oplot, x,y
	  x = x / 2
	  ;oplot, x,y
  end ; if PLOT
end  ; if FIT 

DONE:
end

;===================================================================================================

pro plt_cl, l, cl, oplot=oplot, xtype=xtype, ytype=ytype, xrange=xrange, yrange=yrange, title=title, thick=thick, line=line
; xtype and ytype for log representation

if not keyword_set(oplot) then begin
  plot, l, (l*(l+1)*cl)/(2*!pi), xrange=xrange, yrange=yrange, xtype=xtype, ytype=ytype, xtitle='l',$
  ytitle='l(l+1).Cl/2!7p!3' , title=title, thick=thick, line=line
endif else  oplot, l, (l*(l+1)*cl)/(2*!pi), thick=thick, line=line
 
end

;===================================================================================================

function im_getcl, imag, polarfft=polarfft, psize_arcmin=psize_arcmin, indl=indl, ind2l=ind2l, winAppo=winAppo, opt=opt

im_isospec, imag, powspec, polarfft=polarfft, psize_arcmin=psize_arcmin, Cl=Cl, indl=indl, ind2l=ind2l, winAppo=winAppo, opt=opt

return, Cl

end

;=====================================

function im_getps, imag, polarfft=polarfft, x=x, winAppo=winAppo, opt=opt

im_isospec, imag, powspec, polarfft=polarfft, x=x, winAppo=winAppo, opt=opt

return, powspec

end

;=====================================

function im_getcmb, Cl=Cl, Npix=Npix, Nsimu=Nsimu, psize_arcmin=psize_arcmin
COMMON C_PLANCK

if not keyword_set(Nsimu) then Nsimu = 1

if not keyword_set(Cl) then BEGIN
  ; t = mrdfits('$MRS/data/wmap_lcdm_bf_model_yr1_v1.fits', 1, Header) 
  ; Cl = t.temperature
   Cl = readfits('$MRS/data/def_cl.fits')
   ; help, Cl
   END
vs = size(Cl)
Ncl = vs[1]
lmax =  Ncl
; if lmax GT P_LMAX then lmax = P_LMAX
; tek_color
; plotcl, Cl

if not keyword_set(Npix) then Npix=512L  
N = long(Npix)

; By default the pixel size is the PLANCK pixel size (nside=2048) in arc minute
if not keyword_set(psize_arcmin) then psize_arcmin=1.71774    ; pixel_size(2048)
print, "Pixel Size (Arcmin) = ", psize_arcmin
print, "Number of pixels     = ", strc(Npix) + "x" + strc(Npix)
print, "Lmax Lambda CDM    = ",  Lmax
x=indgen(N/2 + 1)/ float(N)

;l normalisation
factor1 = convert_l2x(psize_arcmin, N)
indl = x* factor1
ind2l = indl*(indl +1)/(2*!pi)

RE = fltarr(N,N)
IM = fltarr(N,N)
U = complexarr(N,N)
Cpt = 0L
TabCMB = fltarr(N,N, Nsimu)
NormSpec = 1. / convert_ps2cl(psize_arcmin, N) / N^2 * 2.

;definition de l'image isotrope via sa FFT
MaxKL=0.
for s=0L, Nsimu-1 do begin
for I=0L,N-1 do begin
	for J=0L,N-1 do begin
	    ix = float(I-N/2.)  / float(N)
	    iy = float(J-N/2.)  / float(N)
  		K=sqrt( ix^2. + iy^2.);distance au centre
  		KL = long(K*factor1+0.5)
  		if KL GT MaxKL then MaxKL = KL
  		if KL GT Ncl-1 then U[I,J]=complex(0.,0.) $
  		else begin
  		  Cpt = Cpt + 1
   		  H=Cl[KL] * NormSpec  ;spectre de puissance a la distance K
  	      RE=randomn(seed)*sqrt(H/2.0)
  	      IM=randomn(seed)*sqrt(H/2.0)
  		  U[I,J]=complex(RE,IM)
  		end
	endfor
endfor
TabCMB[*,*,s] = float( dfti(u)*N_ELEMENTS(u))
endfor

print, "Max used L     = ",  MaxKL

return, TabCMB
end

 

;=====================================


pro ex_im_cmb
; TEST 1
window, 0
window, 1
psize_arcmin=0
cmb = im_getcmb(Cl=Cl, Npix=N, psize_arcmin=psize_arcmin)
info, cmb
help, N
vs = size(Cl)
Ncl = vs[1]
indl = indgen(ncl)

print, psize_arcmin
wset, 0
clcmb = im_getcl(cmb, indl=l1, psize_arcmin=psize_arcmin, /winAppo)
help, clcmb
plot, l1, clcmb*l1*(l1+1)/(2.*!pi), xrange=[0,3000]
oplot, indl, Cl*indl*(indl+1) /(2.*!pi)
wset, 1
plot, indl, Cl*indl*(indl+1) /(2.*!pi)

save, filename='XXCL.xdr', Cl, cmb, clcmb, indl, l1
end

;=====================================

pro ttcl2

; TEST 2
nside=512
N=512.
psize_arcmin = pixel_size(512)
s = mrs_getcmb(cl=cl, nside=512)
F = H2F(s)
F1 = F[*,*,0]
ClF1 = im_getcl(F1, indl=lf1, psize_arcmin=psize_arcmin, /polar)
plot, lf1, ClF1* lf1*(lf1 +1)/(2.*!pi), xrange=[0,2000]
oplotcl, cl
end

;=====================================

