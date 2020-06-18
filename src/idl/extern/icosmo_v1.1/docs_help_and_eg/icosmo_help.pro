pro icosmo_help,routine,list=list

; Aug 08 - modified by AA to read header of the .pro files
; Written by Adam Amara 29th May 2008
;
; PURPOSE: This routine displays help text for the main icosmo
; routines
; INPUT: routine: name of the routine (string)
; OPTIONAL INPUT: list: this will list all the routines that have an
; icosmo_help entry.
; EXAMPLE: icosmo_help,'mk_cosmo'
; -----
;

; -- Construct a list of the routines that have a help entry --
list_routines=['set_fiducial','mk_cosmo','mk_survey','mk_cl_tomo','mk_fisher_lens','mk_bao','plt_cl','prt_precision','asinh','mk_sh']



if (not keyword_set(list)) and (not keyword_set(routine)) then begin
   print,'########################################################################'
   print,'iCosmo_help is a routine that prints help advice for the routines in '
   print,'the iCosmo package. For a list of routines covered in this help file '
   print,'type: icosmo_help,/list'
   print,''
   print,'For help on a specific routine (for example mk_evol) use the following'
   print,'syntax: icosmo_help,''set_fiducial'' ' 
   print,'########################################################################'
   return
endif 

; If the keyword list is set then icosmo_help list all the routines
; that are covered in the help file.
if keyword_set(list) then begin
   print,'########################################################################'
   print,''
   print,'iCosmo_help has help entries for the following routines:'
   print,'---------------------------------------------------------'
   print,transpose(list_routines)
   print,''
   print,'########################################################################'
endif

if not keyword_set(routine) then return

; -- Check that the input routine string does not end in .pro, if so
;    .pro is removed --
if (strmid(routine,3,/reverse) eq '.pro') then routine=strmid(routine,0,strlen(routine)-4)

temp=''
flag=0
help_text='Sorry, This routine does not have a help entry.'
file=file_which(routine+'.pro',/INCLUDE_CURRENT_DIR)
if not keyword_set(file) then begin
   print,'; ***************// HELP BEGIN //**************'
   print,'The routine '+routine+' could not be found.'
   print,'; ***************// HELP END //**************'
   return
endif

openr,10,file
;      while ~ EOF(10) do begin
while not (strmatch(temp,'*// HELP BEGIN //*') or EOF(10)) do begin
   readf,10,temp
endwhile
print,'; ***************// HELP BEGIN //**************'
while not (strmatch(temp,'*// HELP END //*') or EOF(10)) do begin
   flag=1
   readf,10,temp
   print,temp
endwhile
close,10
if (flag eq 0) then begin
   print,'This routine does not have a help entry. Source file is:'
   print,file
   print,'; ***************// HELP END //**************'
endif
end



