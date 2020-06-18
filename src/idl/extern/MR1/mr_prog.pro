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

pro MR_PROG, PROG, DataIn, DataOut,  OPT=Opt, Sphere=Sphere, Verb=Verb
COMMON C_PLANCK

DataIn = float(DataIn)
vsize = size(DataIn)
if N_PARAMS() LT 3 then begin 
        print, 'CALLING SEQUENCE: MR_PROG, PROG, DataIn, DataOut, OPT=Opt'
        goto, DONE
        end

if vsize(0) EQ 0 then begin
        print, 'ERROR: bad first parameter ...'
        print, 'CALLING SEQUENCE: IM_PROG, PROG, DataIn, DataOut, OPT=Opt'
         goto, DONE
        end
 
if not keyword_set(Opt) then Opt = ' '

filename = gettmpfilename()   ; strcompress('xx_temp_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)

; p = strpos(filename, '.fits')
; if p LT 0 then filename = filename + '.fits'

NameData  = gettmpfilename()   ; strcompress('xx_tempdata_'+string(floor(1000000*randomu(seed)))+'.fits',/remove_all)

if not keyword_set(Sphere) then writefits, NameData, DataIn $
else mrs_write, NameData, DataIn

; print, "WRITE ", NameData, " " , P_WAIT
wait, P_WAIT

com = PROG + ' ' + Opt + ' ' + NameData + ' ' +  filename
if keyword_set(Verb) then print, com
spawn, com
; help, filename
; print, "READ"
if not keyword_set(Sphere) then  DataOut = readfits(filename) $
else DataOut =  mrs_read(filename)
; help, dataout

delete, NameData
delete, filename
DONE:

END
