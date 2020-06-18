pro mrs_betatrans, dataIn, Out, wtrans, NbrScale=NbrScale, scaling=scaling, exp=exp

if  not keyword_set(NbrScale) then NbrScale = 6
if NbrScale le 1 or NbrScale ge 20 then begin print,'Error: Number of scales should be between 2 and 20'
    	    	    	    	    goto, DONE
				    end

if keyword_set(exp) then data = exp(datain) else data=datain				    
npix = (size(Data))[1]
nside = npix2nside(npix)
TabWavelet = fltarr(Npix, NbrScale)				    

HScale = Data
ns = nside			   
for j=0,NbrScale-2 do begin
  ns1 = ns / 2
  LScale = mrs_resize(Hscale, nside=ns1)
  LHScale = mrs_resize(Lscale, nside=nside)
  if not keyword_set(scaling) then TabWavelet[*,j] = HScale - LHScale $
  else TabWavelet[*,j] = HScale
  Hscale = LHScale
  ns = ns / 2
end

TabWavelet[*,NbrScale-1] =  LHScale

TabNorm=[0.85,0.12,0.046,0.0224606,0.011,0.006] 
out = {UseGLESP: 0, NbrScale : NbrScale, nside : nside, nx: 0, np:0, npix:npix, Coef : TabWavelet, lmax:0, MeyerWave:0, $
       DifInSH:0, pyrtrans:0, x_sky:0,  y_sky :0, TabNorm:TabNorm, Healpix_with_Glesp: 0}
DONE:


END
