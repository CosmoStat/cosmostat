;+
; NAME:
;        df_get_redshift
;;
; PURPOSE:
;       Calculate the redshifts of a spectra or an array of spectra, using an array of templates. For a large number of spectrum, it is better to call the routine
;       
;
; CALLING:
;     Redshift =  df_get_redshift(Spectrum_or_TabSpectrum, TabTemplate,  lstep= lstep, shiftpar=shiftpar, InfoTemp=InfoTemp)
;
; INPUTS:
;     Spectrum_or_TabSpectrum -- IDL 1D or 2D array :  input Spectra or  TabSpectra[*, 0:N]   N input spectra
;     TabTemplate-- IDL  2D array :   inpute template TabTemplate[*,0:T] 
;    
; OUTPUTS:
;     Redshift -- double or  IDL 1D array :  redshift of the N spectra
;
; INPUT KEYWORDS:
;  lstep -- double : this value is the width of a single pixel on the log-lambda axis; SDSS call this 'Log10 dispersion per pixel'
;  shiftpar -- int: difference in pixels between the end of the template and the end the input spectrum
;
; INPUT/OUTPUT KEYWORDS:
 ;   InfoTemp -- IDL structure: if the routine is called many thimes with the same templates and the same kind of spectrum (i.e. wavelength range, number of pixels, etc),
 ;                        then temporaty array are stored in InfoTemp that can be reused for the next call.
;
; EXAMPLE:
;       Calculate the redshidt of a of spectra
;               TabRedshift =  df_get_redshift(TabSpectra, TabTemplate)
;         
; HISTORY:
;	Written: Daniel Machado & Jean-Luc Starck, 2013
;	Sept, 2013 File creation

;--------------------------------------------------------------------------------------------------------
 
;===========================================

function get_padfft, S,  conj=conj, TargetN
vs = size(S)
N = vs[1]
if TargetN LT N then TargetN = N
DN = (TargetN - N) / 2
Sp = dblarr(TargetN)
Sp[0:N-1] = S
Sp = Shift(Sp, DN)
Spc = dcomplexarr(TargetN)
Spc[*] = Sp[*]
FT = dft(Spc)
if keyword_set(conj) then FT = conj(FT)
return, FT
end

;===========================================

function get_conv_padfft, S, TF
vs = size(TF)
Nt = vs[1]
C = get_padfft(S, Nt)
C1 = C * TF

return, real_part(dfti(C1))

end

;===========================================

function df_get_redshift_one_spectrum, g,  Template, lstep= lstep, chisquare=chisquare, shiftpar=shiftpar, InfoTemp=InfoTemp, MaxChiSquare= MaxChiSquare
; g = denoissed data, clean of baseline
; Template = eigentemplate 


if N_PARAMS() LT 2  then begin 
        lzest =-1
        print, 'CALLING SEQUENCE: Redshift =  df_get_redshift(Spectrum, TabTemplate,  lstep= lstep, shiftpar=shiftpar, InfoTemp=InfoTemp)'
        goto, DONE
        end

if not keyword_set (lstep) then lstep = 1.0E-4 ;this value is the width of a single pixel on the log-lambda axis; SDSS call this 'Log10 dispersion per pixel'
if not keyword_set(shiftpar) then shiftpar = 0
;cross-correlation maxima IT IS BEST TO CHECK THIS VALUE. In general you want your eigentemplates to be longer than your spectra.

vs = size(g)

; if the a full catalog is given, call another routine to estimate all redshift in once
if vs[0] GT 1 then begin
       lzest =-1
        print, 'Error: input spectrum array does not have a correct dimension (1D)… '
        print, 'CALLING SEQUENCE: Redshift =  df_get_redshift_one_spectrum(Spectrum, TabTemplate,  lstep= lstep, shiftpar=shiftpar, InfoTemp=InfoTemp)'
        goto, DONE
 end
      
; Now we have here only one spectrum
 
NbrTemplate = (size(Template))[2]  
NxSpectrum =  (size(g))[1] 
NxTemplate = (size(Template))[1]

ZeroPad_Spectrum = [replicate(0,2* NxTemplate + NxSpectrum), g] ; zero-padding of spectra for cross correlation
NxPadSpectrum = (size(ZeroPad_Spectrum))[1]
FT_Spectrum = fft(ZeroPad_Spectrum, dimension=1,/double) ;fourier transform of treated data
if not keyword_set(InfoTemp) then begin
	ZeroPad_Template = [replicate(0, NxTemplate +2* NxSpectrum, NbrTemplate),   Template] 
;plot, zeropad_template[*,0]
;stop; zero-padding of eigentemplates for cross-correlation
	NxPadTemplate = (size(ZeroPad_Template))[1]
	NxPadMin = min([NxPadSpectrum, NxPadTemplate])
	Conj_FT_Template = conj(fft(ZeroPad_Template, dimension=1,/double)) ;complex conjugate of the fourier transform of the eigentemplates
	Re_Conj_FT_Template = real_part(Conj_FT_Template) ;separating out the real and imaginary parts.
	Im_Conj_FT_Template = imaginary(Conj_FT_Template)
	InfoTemp = {NbrTemplate : NbrTemplate , Re_Conj_FT_Template: Re_Conj_FT_Template,  Im_Conj_FT_Template: Im_Conj_FT_Template, NxPadMin: NxPadMin, NxPadTemplate: NxPadTemplate}
end else begin
	NxPadMin = InfoTemp. NxPadMin
	NxPadTemplate = InfoTemp. NxPadTemplate
	Re_Conj_FT_Template = InfoTemp. Re_Conj_FT_Template
	Im_Conj_FT_Template = InfoTemp. Im_Conj_FT_Template
end
;plot,zeropad_template

prechiterms = dblarr(NxPadMin, NbrTemplate, /nozero); make a large array to hold all the terms contributing to the chisquare per galaxy spectrum
for p = 0, NbrTemplate-1 do begin  
      ; bjzsq=(shift((real_part( fft((gnk * (dcomplex(rebin(rcejk[*,p],padbinsize[0],numbg),rebin(icejk[*,p],padbinsize[0],numbg)))),dimension=1,/inverse, /double)))^2,(-shiftpar),0))
      ReRebin = rebin(Re_Conj_FT_Template[*,p],   NxPadTemplate)
     ;  if p EQ 0 then info, real_part( ReRebin)
      ImRebin =  rebin(Im_Conj_FT_Template[*,p],  NxPadTemplate)
      CF_Rebin = dcomplex(ReRebin,   ImRebin ) * FT_Spectrum
      Cor = real_part( fft(CF_Rebin,  dimension=1, /inverse, /double))
      Cor2 = Cor^2
      bjzsq=shift(Cor2, -shiftpar)
     prechiterms[*,p] = bjzsq 
endfor
if nbrtemplate eq 1 then chisquare = prechiterms else $ 
   chisquare = total(prechiterms,2)
;print,size(chisquare)
;read,idum

; chisquare = dblarr(NxPadMin)
; for i=0, NxPadMin-1 do  chisquare(i) = max(prechiterms[i,*])
lzest =0
lshiftlist = 0
 
 ;maximum of the chisquare.
 MaxChiSquare = max(chisquare[0:NxTemplate-1])
lshiftlist = double(where(chisquare[0:NxTemplate-1] eq MaxChiSquare, c))  ;find the bin number corresponding to the

lzest = double((10d ^ (lshiftlist * Lstep)) - 1d ) ; calculate the redshift estimates from the bin shifts.
chisquare=chisquare[0:NxTemplate-1] 
;info,chisquare
;read,idum

DONE:
  return, lzest  
  
END
 

 ;===========================================

function df_get_redshift, g,  Template, lstep= lstep, chisquare=chisquare, shiftpar=shiftpar, InfoTemp=InfoTemp, TabMaxChiSquare = TabMaxChiSquare
; g = denoissed data, clean of baseline
; Template = eigentemplate 

if N_PARAMS() LT 2  then begin 
        TabZ =-1
        print, 'CALLING SEQUENCE: Redshift =  df_get_redshift(Spectrum, TabTemplate,  lstep= lstep, shiftpar=shiftpar, InfoTemp=InfoTemp)'
        goto, DONE
        end

if not keyword_set (lstep) then lstep = 1.0E-4 ;this value is the width of a single pixel on the log-lambda axis; SDSS call this 'Log10 dispersion per pixel'
if not keyword_set(shiftpar) then shiftpar = 0

vs = size(g)
if vs[0] LT  1 or vs[0] GT  2 then begin
        TabZ =-1
        help, g
       print, 'Error: input spectrum array does not have a correct dimension ... '
        goto, DONE
end
vt = size(Template)
if  vt[0] NE  2 then begin
        TabZ =-1
        help, Template
       print, 'Error: input template array does not have a correct dimension ... '
        goto, DONE
end


if vs[0] EQ 1 then NbrGal = 1 $
else NbrGal =  vs[2] 

TabZ = dblarr(NbrGal)
InfoTemp=0
chisquare=0
TabMaxChiSquare= dblarr(NbrGal)
for i=0, NbrGal-1 do begin
   if total(abs(g[*,i])) gt 0 then begin
      TabZ[i]  = df_get_redshift_one_spectrum(reform(g[*,i]),  Template, lstep= lstep, shiftpar=shiftpar, InfoTemp=InfoTemp, MaxChiSquare= MaxChiSquare, chisquare=chisquare)
      TabMaxChiSquare[i] = MaxChiSquare
   endif else TabZ[i] = -1
 end
DONE:
Infotemp=chisquare
  return, TabZ  
  
END
 
  ;===========================================

