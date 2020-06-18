;-------------------------------------------------------------------------------
;
; NAME:
;        ATWT2D
;
; PURPOSE:
;       Computes the a trous wavelet transform of an image. 
;       The output is a 3D IDL array. 
;
; CALLING:
;
;      ATWT2D, Imag, DataTransf, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord
;       
; INPUTS:
;     Imag -- 2D IDL array: image we want to transform
;
; OUTPUTS:
;     DataTransf -- 3D IDL array: Wavelet Transform
;                                 DataTransf(*,*,i) = ith band of the 
;                                                    wavelet transform
;  
; KEYWORDS:
;      Nscale -- int: Number of scales. Default is 4.
;                     DataTransf is a cube of Nscale planes: 
;                                      DataTransf(*,*,0: Nscale-1)
;			There is no test of the validity of the number of scales with respect
;			to the size of the input images. 
;
;      TabStat -- 2D IDL Array: Statistic of the transform
;                               TabStat[i,0] = mean value of band i
;                               TabStat[i,1] = sigma value of band i
;                               TabStat[i,2] = skewness value of band i
;                               TabStat[i,3] = kurtosis value of band i
;
;      NoBord -- int: if NoBord is set, the wavelet coefficients 
;                     which present a border problem are not taken into
;                     account in the statistics calculation.
; 
; EXAMPLE:
;
;       Compute the multiresolution of the image I with default options
;       (i.e. a trou algorithm with 4 scales).  
;               ATWT2D, I, Output
;
; HISTORY:
;       Written: Jean-Luc Starck 2002.
;       June, 2002 File creation
;
; REFERENCE:
;    J.L. Starck and F. Murtagh, 
;    "Image Restoration with Noise Suppression 
;    Using the Wavelet Transform",
;    Astronomy and Astrophysics, 288, pp-343-348, 1994.
;
; AUTHOR:
;    J.L. Starck
;    Service d'Astrophysique, Centre d'Etudes de SACLAY,
;    Orme des Merisiers, 91191 GIF-Sur-YVETTE CEDEX, France 
;    Email: jstarck@cea.fr        Tel:  +33 (0) 169 085 764
;    http://jstarck.free.fr       Fax:  +33 (0) 169 086 577
;-------------------------------------------------------------------------------
;-------------------------------------------------------------------------------

function test_ind, ind, N

	;this function implements mirror like limit conditions on the edges of the input image
	;ATTENTION : the output may still be out of range ie not in [0, N-1]
	;refinements may be necessary, although not meaningful in practice
	
	ret = ind
	if ind LT 0 then ret = -ind $
	else if ind GE N then ret = 2*N-ind-2
	
	return, ret
	
end
;-------------------------------------------------------------------------------

pro b3splineFast, Im_in, Im_Out,  Step
; /edge_wra
   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   KSize = 4*Step+1
   KS2 = KSize/2
   Kernel = fltarr(KSize)
   Kernel[0] = C1
   Kernel[KSize-1] = C1
   Kernel[KS2+Step] = C2
   Kernel[KS2-Step] = C2
   Kernel[KS2]=C3
   z = convol(Im_in,Kernel, /EDGE_TRUNCATE) 
   kernelY = transpose(Kernel)
   Im_Out = convol(z, kernelY, /EDGE_TRUNCATE)
end

;-------------------------------------------------------------------------------

pro b3spline, Im_in, Im_Out,  Step
   
   Nl = (size(Im_in))(1)
   Nc = (size(Im_in))(2)

   C1 = 1./16.
   C2 = 1./4.
   C3 = 3./8.
   Buffs = fltarr(Nl,Nc)

   for i = 0, Nl-1 do begin
    for j = 0, Nc-1 do begin
      jm = test_ind(j-Step, Nc)
      jp = test_ind(j+Step, Nc)
      jm2 = test_ind(j-2*Step, Nc)
      jp2 = test_ind(j+2*Step, Nc)
      Buffs (i,j) =$ 
		C3 * Im_in(i,j) + C2 * (Im_in(i,jm) + Im_in(i,jp))$
		+ C1 * (Im_in(i,jm2) + Im_in(i,jp2)) 
	endfor
   endfor

   for i = 0, Nl-1 do begin
    for j = 0, Nc-1 do begin
      im = test_ind(i-Step, Nl)
      ip = test_ind(i+Step, Nl)
      im2 = test_ind(i-2*Step, Nl)
      ip2 = test_ind(i+2*Step, Nl)
      Im_Out (i,j) =$
			C3 * Buffs(i,j) + C2 * (Buffs(im,j) + Buffs(ip,j))$
			+ C1 * (Buffs(im2,j) + Buffs(ip2,j))
    endfor
   endfor
end

;-------------------------------------------------------------------------------


pro atwt2d, Ima, WT, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord, modif=modif, nostat=nostat

	;-------------------------------------------------------------------------------
	;a few verifications
	if N_PARAMS() LT 2 then begin 
        print, 'CALLING SEQUENCE: atwt2d, Ima, WT, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord'
        goto, DONE
	end

	if not keyword_set(Nscale) then Nscale =4

	if Nscale LT 2 then  Nscale =4
	
	vsize = size(Ima)
	if vsize[0] NE 2 then begin
        print, 'Error: First parameter is not an image ...'
        print, 'CALLING SEQUENCE: atwt2d, Ima, WT, Nscale=Nscale, TabStat=TabStat, NoBord=NoBord'
        goto, DONE
	end
	;-------------------------------------------------------------------------------
	;recovering a few parameter values, initializing the transform loop ...
	Nx = vsize[1]
	Ny = vsize[2]
	Nz = Nscale
	NStep=Nscale-1
	Nstat = 4
	TabStat= fltarr(Nstat,Nscale)
	WT = fltarr(Nx,Ny,Nz)
	Im_in = Ima
	Im_Out = fltarr(Nx,Ny)
	Step_trou = 1

	for i=0, NStep-1 do begin
		;B3 Spline smoothing
		b3splinefast, Im_in, Im_Out, Step_trou 
                if keyword_set(modif) then begin
		   Im_Aux = Im_Out
		   b3splinefast, Im_Out, Im_Aux, Step_trou 
		   WT(*,*,i) = Im_in - Im_Aux
		end else begin
		; Wavelet coefficient calculations
		  WT(*,*,i) = Im_in - Im_Out
		end
		Im_in = Im_Out

		; Statistic calculation
		if not keyword_set(nostat) then begin
		  BordSize = 2*Step_trou
		  if not keyword_set(NoBord) then BordSize = 0
	
		  band = reform(WT(BordSize:Nx-1-BordSize,BordSize:Ny-1-BordSize,i))
		  TabStat[0,i] = total(band, /double) / N_ELEMENTS(band)
		  TabStat[1,i] = sigma(band)  
		  TabStat[2,i] = skewness(band) 
		  TabStat[3,i] = kurtosis(band) 
		  ;print, "Scale ", i+1, "   Sigma = ", TabStat[1,i], "   Skewness = ", TabStat[2,i], "   Kurtosis = ", TabStat[3,i]
                end
		
		; New distance between two pixels
		if i NE NStep-1 then Step_trou = Step_trou * 2
        endfor


	; Smooth array
	i =  NStep  
	WT(*,*,i) = Im_Out
	
	if not keyword_set(nostat) then begin
	  band = reform(WT(BordSize:Nx-1-BordSize,BordSize:Ny-1-BordSize,i))
	  TabStat[0,i] = total(band, /double) / N_ELEMENTS(band)
	  TabStat[1,i] = sigma(band)  
	  TabStat[2,i] = skewness(band) 
	  TabStat[3,i] = kurtosis(band) 
	end
	;print, "Smooth Ima:      Sigma = ", TabStat[1,i], "   Skewness = ", TabStat[2,i], "   Kurtosis = ", TabStat[3,i]

DONE:

end

;-------------------------------------------------------------------------------

pro atrec2d,  WT, RecIma, Modif=Modif
vs = size(WT)
Nscale = vs[3]
NStep=Nscale-1
Step_trou = 1
for i=0, NStep-2 do Step_trou = Step_trou * 2
	
RecIma = WT[*,*,Nscale-1]
if keyword_set(modif) then Im_Out = fltarr(vs[1],vs[2])
for j=Nscale-2,0,-1 do begin
  ; print, j, Step_trou
  if keyword_set(modif) then  begin
      b3splinefast, RecIma, Im_Out, Step_trou
      RecIma = Im_Out + WT[*,*,j]
  end else RecIma = RecIma + WT[*,*,j]
  Step_trou = Step_trou / 2
end
end


;-------------------------------------------------------------------------------




