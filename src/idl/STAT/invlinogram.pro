;+
; NAME:
;	INVLINOGRAM
;	
; PURPOSE:
;	transformee de radon inverse moyennee
;
; CALLING SEQUENCE:
;	INVLINOGRAM, trans, recons, tab_indices = tab_indices
;
; INPUTS:
;	trans -- 2D IDL array: input Radon transform	
;
; INPUT KEYWORDS:
;		tab_indices -- array : if set, this is an index array to be used for the line extraction in the linogram 
;							if not set, then the array is computed using lino_grid.pro
;
; OUTPUTS:
;	recons -- 2D IDL array: output reconstructed image
;
; EXAMPLE:
;       linogram, data, trans
;	invlinogram, trans, recons
;
; HISTORY:
;	Written:Ludovic POUPARD  & Yassir MOUDDEN, 2005.
;       May, 2005 File creation
;------------------------------------------------------------------------


PRO invlinogram,  trans, recons, tab_indices = tab_indices

;------------------------------------------------------------------------
; declarations
;------------------------------------------------------------------------
Ntx = (size(trans))[1]
Nty = (size(trans))[2]
N=Nty
cc=(N-1)/2

ligne_rec=dcomplexarr(2*N-2,N)

TF_image_rec=dcomplexarr(N,N)
recons=dblarr(N,N)
datacount=dblarr(N,N)



;------------------------------------------------------------------------
;reconstruction
;------------------------------------------------------------------------


ligne_rec= fft(trans, -1, dimension = 2) * N_ELEMENTS(trans(0,*))

if not keyword_set(tab_indices)$
then begin
	lino_grid, N, tab_indices
endif

for num_x = 0, Ntx-1 do begin
	temp = tab_indices(num_x, *)
	TF_image_rec(temp) = TF_image_rec(temp) +  ligne_rec(num_x, *)
	datacount(temp)=datacount(temp)+1	
endfor

datacount=datacount+(datacount eq 0)
TF_image_rec=TF_image_rec/datacount


;---------------------------------------------------------
;image reconstruite
;----------------------------------------------------------

recons=REAL_PART(  fft( TF_image_rec,  +1, /double  ) )/N_ELEMENTS(TF_image_rec)	

END

