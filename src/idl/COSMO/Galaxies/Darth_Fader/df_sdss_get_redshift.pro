; zans = zcompute(objflux, objivar, starflux, [starmask, $
;    fixed_template=, nfind=, $
;    poffset=, pspace=, pmin=, pmax=, mindof=, width=, minsep=, $
;    plottitle=, /doplot, /debug, /verbose, zans_fixed= ]
;  poffset    - Offset between all objects and templates, in pixels.
;                A value of 10 indicates that STARFLUX begins ten pixels
;                after OBJFLUX, i.e. OBJFLUX[i+10] = STARFLUX[i] for the
;                case when the relative redshift should be zero.  If the
;                wavelength coverage of the templates is larger, then the
;                value of ZOFFSET will always be negative.
;                [Scalar or vector of size NOBJ]



function df_sdss_get_redshift, g,  Template, lstep= lstep,  diffpar=diffpar, rms=rms, tred=tred
; g = denoissed data, clean of baseline
; Template = eigentemplate 

if N_PARAMS() LT 2  then begin 
        TabZ =-1
        print, 'CALLING SEQUENCE: Redshift =  df_get_redshift(Spectrum, TabTemplate,  lstep= lstep, diffpar=diffpar)'
        goto, DONE
        end

if not keyword_set (lstep) then lstep = 1.0E-4 ;this value is the width of a single pixel on the log-lambda axis; SDSS call this 'Log10 dispersion per pixel'
if not keyword_set(diffpar) then diffpar = 0

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
training_lmin = 3000.  ; Amstrom units
data_lmin = 3600.
poffset =  round((alog10(training_lmin) - alog10(data_lmin))/lstep)

for i=0, NbrGal-1 do begin
    TabZ[i]  = df_get_redshift_one_spectrum(g[*,i],  Template, lstep= lstep, diffpar=diffpar)
    if keyword_set(rms) then objivar = 1. / rms^2 $
    else objivar = g[*,i] * 0 + 1.
    zans = zcompute(g[*,i], objivar, Template, poffset=poffset)
    hs, zans
    if keyword_set(Tred) then print, "True-z = ", Tred, ", Est-z", TabZ[i] $
    else print,  TabZ[i]
    lshiftlist = zans.z
    lzest = double((10d ^ (lshiftlist * Lstep)) - 1d ) 
    help, lzest
 end
DONE:
  return, TabZ  
  
END
 
  ;===========================================

