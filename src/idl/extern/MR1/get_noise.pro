FUNCTION GET_NOISE, Data, Niter=Niter
;+ 
; NAME: 
;     GET_NOISE
;
; PURPOSE: 
;    Find the standard deviation of a white gaussian noise in the data.
;
; CALLING SEQUENCE: 
;   output=GET_NOISE(Data)
;
; INPUTS: 
;   Data -- IDL array
;
; OPTIONAL INPUT PARAMETERS: 
;   none
;
; KEYED INPUTS: 
;   Niter --scalar: number of iterations for k-sigma clipping
;                   default is 3.
;
; OUTPUTS: 
;    output
;
; MODIFICATION HISTORY: 
;    17-Jan-1996 JL Starck written with template_gen 
;-
 
;------------------------------------------------------------
; parameters check
;------------------------------------------------------------
 
 IF N_PARAMS() LT 1 THEN BEGIN
   PRINT, 'CALLING SEQUENCE: ', 'output=GET_NOISE(Data)'
   GOTO, CLOSING
 ENDIF
 
;------------------------------------------------------------
; function body
;------------------------------------------------------------

if not keyword_set(Niter) then Niter = 3.
vsize = size(Data)
dim = vsize[0]
sigma=-1
case dim of 
  3: BEGIN
        nco=vsize[1]
        nli=vsize[2]
        npz=vsize[3]
        indices=indgen(npz-2)+1
        D_cube=fltarr(nco,nli,npz)
        c1=-1./sqrt(6.)
        c2=2./sqrt(6.)
        D_cube[*,*,1:npz-2]=c1*(Data[*,*,indices-1] + Data[*,*,indices+1]) + $
		            c2*(Data[*,*,indices])
        D_cube[*,*,0]=c2*(Data[*,*,0]-Data[*,*,1])
        D_cube[*,*,npz-1]=c2*(Data[*,*,npz-1]-Data[*,*,npz-2])
        sigma=sigma_clip(D_cube,Niter=Niter)
     END
  2: BEGIN
        im_smooth, Data, ima_med, winsize=3, method='median'
        sigma = sigma_clip(Data-ima_med,Niter=Niter) / 0.969684
     END
  1: BEGIN
        Sigma = sigma_clip(Data - median(Data,3), Niter=Niter) / 0.893421
     END
  else: BEGIN
        print, 'Error: bad dimension of the input parameter'
        END
 END

;------------------------------------------------------------
; closing
;------------------------------------------------------------
 
 CLOSING:
 
  RETURN, Sigma
 
 END
