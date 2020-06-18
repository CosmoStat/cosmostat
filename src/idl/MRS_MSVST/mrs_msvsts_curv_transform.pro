;+
; NAME:
;        MRS_MSVSTS_CURV_TRANSFORM
;
; PURPOSE:
;	    Compute the multi-scale variance stabilising transform
;       on the sphere with standard undecimated curvelet transform on the sphere, using the healPix pixel representation 
;       (nested data representation). A band of the curvelet transform is defined by two number, the 
;       2D WT scale number and the ridgelet scale number.
;       The output is a IDL structure.
;       A band at wavelet scale j (j=0..NBRSCALE-1) and ridglet scale j1 can be extracted using the  
;       function mrs_curget(Curtrans, j, j1) (ex: Scale2_1 = mrs_curget(CurTrans, 2, 1))
;       and a band can be inserted in the transformation using the routine  mrs_curput
;       (ex:  mrs_curput, CurTrans, Scale2_1, 2, 1).
;     
;
; CALLING:
;
;     mrs_msvsts_curv_transform, Imag, Trans, lmax=lmax, NbrScale=NbrScale, FirstBlockSize=FirstBlockSize
;
; INPUTS:
;     Imag -- IDL array of healpix map: Input image be transformed 
;    
; OUTPUTS:
;     Trans -- IDL structures with the following fields:  
;    NBRSCALE      --   INT: Nbr of the scale in the 2D WT
;    TABBLOCKSIZE  --   INT: TABBLOCKSIZE[j], Block size in the ridgelet transform at scale j  
;                                             j = [0..NBRSCALE-2]
;    TABNBRSCALERID --  INT:  TABNBRSCALERID[j], number of ridgelet band at scale j 
;    TABNORM       --   2D IDL ARRAY: Normalization array
;    RIDSCALE1     --   IDL STRUCT: ridgelet transform of the first wavelet scale
;                                   (see mrs_ridtrans.pro for details)
;      ...
;    RIDSCALEj     --   IDL STRUCT : ridgelet transform of the jth wavelet scale
;                                    j = [0..NBRSCALE-2]
;    LASTSCALE     --   IDL 1D array: Healpix image of the coarsest scale
;    WT            --   IDL STRUCT: Wavelet structure (for internal use only)
;    PYRTRANS      --   INT: equal to 1 for a pyramidal curvelet transform and 0 otherwise 
;
; KEYWORDS:
;      NbrScale -- INT: Number of scale in the 2D wavelet transform (defaut 4)
;      Undec -- INT: if set, an undecimated curvelet transform is used instead of the pyramidal
;                    curvelet transform
;      FirstBlockSize -- INT: Block size in the ridgelet transform at 
;                             the finest scale (default is 16)
;      Lmax      : Number of used spherical harmoniques used in the wavelet transform
;                 (defaut = 3*nside, should be between 2*nside and 4*nside)
;      Overlap   -- LONG: is equal to 1 if blocks are overlapping
;      
;
; EXTERNAL CALLS:
;        mrs_ridtrans (IDL) program
;
; EXAMPLE:
;
;       Compute the curvelet transform of an image I with default options
;        The result is stored in Output
;               mrs_msvsts_curv_transform, Imag, Output 
;         
; HISTORY:
;	Written:  Jérémy Schmitt, 2010
;	June, 2010 File creation
;-----------------------------------------------------------------


;-----------------------------------------------------------------

pro mrs_cursig, Trans
for s2d = 0,Trans.NbrScale-2 do begin
  for s1d = 0,Trans.TabNbrScaleRid[s2d]-1 do begin
    Band =  mrs_curget(Trans, s2d, s1d) / sqrt( Trans.TabBlockSize[s2d])
    print, s2d, s1d, sigma(Band)
  end
end

end
;-----------------------------------------------------------------

pro mrs_msvsts_curv_transform, Imag, CurTrans, Opt=Opt, lmax=lmax, NbrScale=NbrScale, Overlap=Overlap, SSR=SSR, FirstBlockSize=FirstBlockSize, Silent=Silent
COMMON MR1ENV

if N_PARAMS() LT 2  then begin 
        print, 'CALLING SEQUENCE: mrs_wttrans, Imag, RidTrans,  Opt=Opt'
        goto, DONE
        end
	 
npix = (size(imag))[1]
nside = npix2nside(npix)

;if keyword_set(undec) then
;wtt_msvst, Imag, WT, NbrScale=NbrScale, lmax=lmax, /Healpix_with_Glesp
mrs_msvsts_iuwt_transform, Imag, WT, NbrScale=NbrScale, lmax=lmax
;else mrs_pwttrans, Imag, WT, NbrScale=NbrScale, lmax=lmax

NS = WT.NbrScale 
TabBlockSize = intarr(NS)
TabNbrScaleRid = intarr(NS)
if  not keyword_set(Opt) then Opt= ' '
if  not keyword_set(FirstBlockSize) then FirstBlockSize = 16
HalfSize = FirstBlockSize / 2
for j=0,NS-1 do begin 
  if keyword_set(SSR) OR NOT keyword_set(MR1OK) then TabBlockSize[j] = 2*HalfSize $
  else TabBlockSize[j] = 2*HalfSize
  if not keyword_set(silent) then print, 'Scale ', j+1, ' BlockSize = ', TabBlockSize[j]
  ;if not keyword_set(undec) then if j mod 2 EQ 1 then HalfSize = HalfSize / 2
  if HalfSize LT 4 then HalfSize = 4
  ;if keyword_set(undec) then 
  if j mod 2 EQ 0 then HalfSize = HalfSize * 2
end

BandName = strarr(NS)
for j=0,NS-2 do begin
    Scale = mrs_wtget(WT, j)
    OptRid = Opt    ;  + ' -b' + strcompress(string(TabBlockSize[j]), /remove_all)
    ; if not keyword_set(silent) then print, 'Ridgelet transform: ', j+1, ' Option = ', OptRid
	mrs_ridtrans, Scale, Rid, Opt=OptRid, BlockSize=TabBlockSize[j], Overlap=Overlap 
	;mrs_ridtrans, Scale, Rid, Opt=OptRid, BlockSize=TabBlockSize[j],Overlap=Overlap,NbrScale=4
	print,'echelles de ridgelet :',Rid.Nbrscale
    TabNbrScaleRid[j] = Rid.NbrScale
    my_command = 'RidScale'+strcompress(string(j+1), /remove_all)
    BandName[j] = my_command 
    my_command = my_command +'= Rid'
    my_command = strcompress( my_command, /remove_all)
    ; print, 'cmd = ',  my_command
    ACK = EXECUTE( my_command)    
end
LastScale =  mrs_wtget(WT, NS-1)

;if keyword_set(undec) then
;TabNorm = [ [0.829881, 0.134254, 0.0948845, 0.0948845, 0.0948845], $
;            [0.0763478, 0.0882548,0.0536306,0.0415498,0.0415498], $
;            [0.0181707, 0.0348571,0.0496051,0.0532982,0.0532982], $
;	    [0.00517699, 0.0102914,0.0226026,0.0592364,0.0592364], $
;            [0.00238540, 0.00377852, 0.00838926, 0.0473952, 0.0473952], $
;            [0.000855462, 0.00125520,0.00264597,0.00612633, 0.0343854]]
;else  TabNorm = [ [0.829880, 0.134252, 0.0948756],$
;                  [0.106301, 0.0434990, 0.0279083],$
;                  [0.0438302, 0.0291678, 0.01],$
;                  [0.0213290,0.0139984, 0.006],$
;                  [0.0110710, 0.00936899, 0.001],$
;                  [0.00592302, 0.00402646, 0.003]]

;TabNorm : tableau des écarts-types de chaque bande de curvelet pour un bruit poissonien d'écart-type 100 (6 echelles d'ondelettes)
;TabNorm = [ [1.74550, 0.348175], $
;            [0.230621, 0.248233, 0.196981], $
;            [0.0548140, 0.0989918, 0.219056], $
;	        [0.0212912, 0.0417454, 0.0875663, 0.20375], $
;            [0.00989616, 0.0158273, 0.0352021, 0.163248], $
;            [0.0777546]]

TabNorm = [ [1.74550, 0.348175,1,1], $
            [0.230621, 0.248233, 0.196981,1], $
            [0.0548140, 0.0989918, 0.219056,1], $
	        [0.0212912, 0.0417454, 0.0875663, 0.20375], $
            [0.00989616, 0.0158273, 0.0352021, 0.163248]]

  pyrtrans = 1
  if keyword_set(undec) then pyrtrans = 0
  my_command = 'CurTrans = { NbrScale : NbrScale,'
  my_command = my_command+'TabBlockSize : TabBlockSize, '
  my_command = my_command+' TabNbrScaleRid : TabNbrScaleRid, '
  my_command = my_command+' TabNorm :TabNorm,'
  for j=0, NS-2 do  my_command = my_command + BandName[j] +':'+ BandName[j] +','
  my_command = my_command + ' LastScale :LastScale, '
  my_command = my_command + ' WT :WT, '
  my_command = my_command+' pyrtrans :pyrtrans'
  my_command = strcompress( my_command, /remove_all)
  my_command =my_command+'}'
  ; print, my_command
  ACK = EXECUTE( my_command)
DONE:

END

