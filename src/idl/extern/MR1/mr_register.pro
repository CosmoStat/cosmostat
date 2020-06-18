;+
; NAME:
;        MR_PROG
;
; PURPOSE:
;	Run a binary program which takes as an input a  fits file
;       and write as an output another fits file.
;
; CALLING:
;      MR_PROG, Prog, DataIn, DataOut, OPT=Opt
;
; INPUTS:
;     Prog -- string: program to be called
;
;     DataIn --  IDL array: input data
;
; OUTPUTS:
;     DataOut --  IDL array: result
;
; KEYWORDS:
;      Opt -- string: string which contains the differents options for the called program.
;
; EXTERNAL CALLS:
;       program called.
;
; EXAMPLE:
;       Compute the 3D PCA transform of a data setr, calling the im3d_pca C++ program.
;               IM3D_PROG, 'im3d_pca', DataIn, PcaOut
;
; HISTORY:
;	Written: Jean-Luc Starck 2004.
;	July, 2004 File creation
;-

pro MR_register, DataIn, DataRef, DataOut,  OPT=Opt

DataIn = float(DataIn)
vsize = size(DataIn)
if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: mr_register, Ima, Ref, RegIma, OPT=Opt'
        goto, DONE
        end

if vsize(0) EQ 0 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: mr_register, Ima, Ref, RegIma, OPT=Opt'
         goto, DONE
        end
 
if not keyword_set(Opt) then Opt = ' '

FN1  = strcompress('xx_temp1_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
FN2  = strcompress('xx_temp2_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)
FN3  = strcompress('xx_temp3_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)

writefits, FN1, DataIn
writefits, FN2, DataRef
com = 'mr_registration ' + Opt + ' ' + FN1 + ' ' +  FN2  + ' ' +  FN3
; print, com
spawn, com
DataOut = readfits(FN3)

delete, FN1
delete, FN2
delete, FN3
DONE:

END
