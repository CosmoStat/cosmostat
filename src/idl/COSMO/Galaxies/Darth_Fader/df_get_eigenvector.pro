;+
; NAME:
;        df_get_eigenvector
;
; PURPOSE:
;   Computes the main eigen vector of a galaxy spectra training set
;
; CALLING:
;     Template =  df_get_eigenvector(Training,  AllEigen=AllEigen,
;     Denoising=Denoising, OptDen=OptDen, EnergPercent=EnergPercent, NTemp=Ntemp)
;
; INPUTS:
;     Training -- IDL 2D array :  Training[*, 0:N-1]   N spectra
;    
; OUTPUTS:
;     Template -- IDL 2D array :  Template[*, 0:T-1]   T spectra
;                 template. Relevant template which contain EnergPercent% of the energy of the data.
;
; INPUT KEYWORDS:
;  Denoising -- Scalar:  if set, a denoising is applied on each spectrum
;  OptFil -- String : option used for the sdenoiing
;  EnergPercent -- double: percentage energy used to define the
;                  relevant eigenvectors. Default is 99.93%.
;  Ntemp -- number of templates to be returned. This overrides
;           EnergPercent and returns the specified number of templates
;
; KEYWORDS:
;      AllEigen   - IDL 2D array :  AllEigen[*, 0:N], all template
;
; EXAMPLE:
;       Compute the spherical harmonix transform of an image. 
;        The result is stored in Output
;               mrs_trans, Imag, Output 
;         
; HISTORY:
;	Written: Daniel Machado & Jean-Luc Starck, 2013
;	Sept, 2013 File creation

;--------------------------------------------------------------------------------------------------------
 
 
;===============================================================

function  df_get_eigenvector, Training,  AllEigen=AllEigen, Denoising=Denoising, OptDen=OptDen, EnergPercent=EnergPercent, NTemp=NTemp

if N_PARAMS() LT 1  then begin 
        erelevant=-1
        print, $
           'CALLING SEQUENCE: Template =  df_get_eigenvector(Training,  AllEigen=AllEigen, Denoising=Denoising, OptDen=OptDen, EnergPercent=EnergPercent, Ntemp=Ntemp)'
        goto, DONE
        end
        
if not keyword_set(EnergPercent) then EnergPercent=99.95d

ntrain = (size(training))[2]  
trains = ntrain-1

tbins = (size(training))[1]

S = Training
;A crude continuum subtraction (at high SNR this procedure is sufficient).

if keyword_set(Denoising) then begin
   if keyword_set(OptDen) then Opt=OptDen $
   else OPT='-M -f3 -n6 -t11 -k -K -i20 -s4 -P'
   mr1d_filter,training,S,OPT=Opt ; these options are explained in mr1d_filter.
;  continuua = training-s
end  

normS = total(S^2d ,1,/double)  ; spectra have to be normalised to the square of the total of each spectrum since they are continuum-subtracted

T = S/(rebin(transpose(sqrt(normS)),tbins,ntrain)) ;divide each spectrum by it's norm.

COMATRIX = T ## transpose(T)         ; create a correlation matrix
A = dblarr(ntrain,ntrain,/nozero)    ; create a blank matrix to fill with eigenvectors.

D = LA_EIGENQL(comatrix,eigenvectors=A,/double) ; routine to diagonalise the correlation matix/calculate eigenvectors and eigenvalues

; Order eigenvalues and eigenvectors in decreasing order of importance:     
R = reverse( transpose(A) )     ; in entries (template) space, this is the right way round to get columns of eigenvectors and
                                ; [transpose(R) ## C ## R] = eigenvalue matrix.

F = reverse(D)                  ; order eigenvalues from largest to smallest.

; Construct the eigentemplates from the eigenvectors:     
E = transpose(R) ## T

percoverage = total(f,/cumul,/double)/total(f,/double)*100d ; calculate the cumulative eigenvalue weight that each eigentemplate carries.

if keyword_set(NTemp) then rel = NTemp $
else rel = min(where(percoverage ge EnergPercent))+1 ; calculate the relevant number of eigentemplates based on how many are required to reach a total
                                ; of 99.93% of the eigenvalue weight.
relm = rel-1
print, "Nbr of templates = ", rel , ", PercentEnergy = ", percoverage[relm]

E = E/(rebin(transpose(sqrt(f)),tbins,ntrain)) ; orthonormalise the eigentemplates such that E ## E = identity.
erelevant = e[*,0:relm]                        ; restrict a new matrix to only the relevant eigentemplates


AllEigen = e

DONE:

return, erelevant

END

