;+
; NAME:
;        mrs_resize
;
; PURPOSE:
;   Resize map either in healpix or glesp representation. If BandLimited is set, a Low-Pass filter is applied to the map. 
;   The low-pass filter is equal to 1 from l=0 to LmaxAnalysis,  and LmaxAnalysis is equal by default to 2*final_nside,
;   and then decrease to 0, at Lmax value, with Lmax = MIN( [3*final_side, 3* LmaxAnalysis]).
;   If a lowpass filtering is given via the keyword  LowPassFilter=LowPassFilter, then the filter is used for the low pass filtering step.
;
; CALLING:
;     NewMap =  mrs_resize,Imag, nside=nside, nx=nx, np=np, ViaAlm=ViaAlm, BandLimited=BandLimited, LmaxAnalysis=LmaxAnalysis, LowPassFilter=LowPassFilter)
;
; INPUTS:
;     Imag -- IDL array of healpix map or GLESP structure: Input image to be transformed 
;    
; OUTPUTS:
;     Trans -- IDL array of healpix map or GLESP structure.
;
; INPUT KEYWORDS:
;      nside     : the nside of the healpix output map
;      nx        : the latitude partition of glesp map
;      np        : the longitude partition of glesp map
;    BandLimited:  if set a low pass filtering is applied
;    LmaxAnalysis: Values until  LmaxAnalysis  are equal to 1 in the low pass filtering, i.e.  LowPassFilter[0: LmaxAnalysis] = 1 
;
; INPUT/OUTPUT KEYWORD:
;     LowPassFilter: fltarr(0:Lmax) -- low pass filtering used in the filtering
;
; EXTERNAL CALLS:
;       anafast (healpix software)
;       cl2map (glesp software)
;       reorder (healpix software)
;
; EXAMPLE:
;       resize an healpix map
;
;               map2 = mrs_resize(map,nside = 256)
;               map3 = mrs_resize(map, nside=16, LmaxAnalysis=5, LowPassFilter=LPout)
;
;       resize an GLESP map
;
;               map2 = mrs_resize(map,nx = 512,np = 1024)
;         
; HISTORY:
;	Written:  Jean-Luc Starck & Pierrick Abrial , 2006
;	February, 2006 File creation
;   Sept 2012, Florent Sureau, add pass band filtering
;--------------------------------------------------------------------------------------------------------


function mrs_resize,Imag, nside=nside, nx=nx, np=np, ViaAlm = ViaAlm, BandLimited=BandLimited, LmaxAnalysis=LmaxAnalysis, LowPassFilter=LowPassFilter
; np = 2nx = lmax/2

    ImagNside = gnside(Imag)
   if not keyword_set(nside) then  nside=ImagNside

   if  not  keyword_set(ViaAlm) then ViaAlm = 0
   if keyword_set(BandLimited) or keyword_set(LmaxAnalysis ) or keyword_set(LowPassFilter )  then ViaAlm = 1
     
    if ImagNside EQ nside then ViaAlm=0
   ; print, ViaALM, ImagNside, nside, BandLimited
   if  ViaAlm EQ 0 then begin
       if type_code(Imag) EQ 8 and keyword_set(nx) and keyword_set(np) then BEGIN
           mrs_almtrans,Imag,alm
           mrs_almrec,alm,map,nx=nx,np=np;,lmax=nx/2
           return,map
      endif else BEGIN 
          if keyword_set(nside) then begin
          ud_grade,Imag,map ,nside_out=nside, order_in='nested',order_out='nested'
          return,map
   endif else return,Imag
   endelse
   return,0
   endif
 
  ; HERE ViaAlm=1. We use the ALM to make the resize
  if not keyword_set(LowPassFilter) then begin
       if keyword_set(LmaxAnalysis ) then BandLimited = 1
      if not keyword_set(LmaxAnalysis) then LmaxAnalysis = 2 * nside - 1
      Lmax = min( [LmaxAnalysis * 3,  3 * nside])
      if keyword_set(BandLimited) then begin
          	beamdata=dblarr(Lmax +1)
	        LowPassFilter =getidealbeam(beamdata, lmin=LmaxAnalysis, lmax=Lmax, /tozero)
     end
  end else begin
     vs = size(LowPassFilter)
     Lmax = vs[1] - 1
  end
   
    NImag = 0.*dblarr(nside^2.*12.,1)
    mrs_almtrans,NImag, A, /tab
    if Lmax  GT a.lmax then Lmax = A.lmax 

    mrs_almtrans,Imag, alm, /tab
    A.alm[*]=0.
    
    if not keyword_set(LowPassFilter) then  A.alm[0:alm.lmax, 0:alm.lmax, *] = alm.alm[0:alm.lmax, 0:alm.lmax, *] $
    else begin
    		for l=0,Lmax  do A.alm[l,0:l,*]=alm.alm[l,0:l,*] * LowPassFilter[l]
    end
    
    mrs_almrec,A, NImag
    return,NImag
        
end
   
    
