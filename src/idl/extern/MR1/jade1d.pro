;-------------------------------------------------------------------------------
;+
; NAME:
;	jade1d
;	
; PURPOSE:
;	Apply the Independant Componant Method to a set of vectors using
;       the JADE method. The routine is identical to jade routine, except
;       that the input-output vectors are ordered differently (V(*,i) for
;       the ith vector, instead of V(i,*)).
;
; EXPLANATION:
;	compute ICA
;
; CALLING SEQUENCE:
;	jade1d, ObservSig, NbSource, DeMixingMat, Process $
;             DispDebug=DispDebug, FileDebug=FileDebug, $
;	      Verbose=Verbose
; INPUTS:
;	ObservSig : mixing of input signal (ObservSig = A # input signals)
;                    ObservSig(*, i) = ith vector
;       NbSource  : number of sources in input signal
;
; OPTIONAL INPUTS:
;	DispDebug : debug trace on screen
;	FileDebug : debug trace on file (Trace.dat in local directory)
;	Verbose	: verbose
;
; KEYWORD PARAMETERS:
;
; OUTPUTS:
;	DeMixingMat : demixing matrix 
;       Process: Reconstructed process = DeMixingMat # ObservSig
;                Process(*, i) = ith Process, with i = 0.. NbSource - 1
;
; DEPENDENCIES:
;
; RESTRICTIONS:
;
; PROCEDURES USED:
;	whitening_signal, estimation_cumulant, init_diagonalisation
;	contrast_optimisation, separating_matrix
;
; EXAMPLE:
;	jade1d, ObservSig, NbSource, DeMixingMat
;
; MODIFICATION HISTORY:
;
;-------------------------------------------------------------------------------   

pro jade1d, ObservSig, NbSource, DeMixingMat, Process, $
          DispDebug=DispDebug, FileDebug=FileDebug, $
	  Verbose=Verbose

if N_PARAMS() LT 3 then begin 
        print, 'CALL SEQUENCE:  jade1d, ObservSig, NbSource, DeMixingMat, Process'
        goto, DONE
        end
	  
Nv = (size(ObservSig))(2)
Np = (size(ObservSig))(1)
tab = dblarr(Nv, Np)
for i=0, Nv-1 do begin
    b = ObservSig(*,i)
    tab(i,*) = b 
    end

jade, tab, NbSource, DeMixingMat, O, $
          DispDebug=DispDebug, FileDebug=FileDebug, $
	  Verbose=Verbose

Process = dblarr(Np, NbSource)
for i=0, NbSource-1 do begin
      b(*) =  O(i,*)
      Process(*,i) = b 
      end
DONE:
end


