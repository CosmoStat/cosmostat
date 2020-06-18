function check_fidinput,cosmo_in=cosmo_in,expt_in=expt_in,calc_in=calc_in

;Routine for checking the user specified input to set_fiducial
; Written by Adam Amara 31 July 2008
; PURPOSE: This routine checks the inputs that the user has specified
; to check that the enteries are valid.
; INPUTS: cosmo_in - the input cosmology options
;         calc_in -  the input calculation options
;         expt_in -  the input experiment options
; OUTPUTS: an integer that is 0 - no error, or 1- error detected.

;create a fiducial structure that is produced with no user
;inputs. This contains all the active tags.
Error_check=0
fid_template=set_fiducial()

if keyword_set(cosmo_in) then begin
   tags=tag_names(cosmo_in)
   for i=0,n_elements(tags)-1 do begin
      if not tag_exist(fid_template.cosmo,tags(i)) then begin
         Error_check=1
         print,'Error: ',tags(i),' is not a valid COSMO input for set_fiducial'
         tags_valid=tag_names(fid_template.cosmo)
         print,'Valid inpts are:'
         print,tags_valid
      endif
   endfor
endif

if keyword_set(calc_in) then begin
  tags=tag_names(calc_in)
   for i=0,n_elements(tags)-1 do begin
      if not tag_exist(fid_template.calc,tags(i)) then begin
         Error_check=1
         print,'Error: ',tags(i),' is not a valid CALC input for set_fiducial'
         tags_valid=tag_names(fid_template.calc)
         print,'Valid inpts are:'
         print,tags_valid
      endif
   endfor
endif

if keyword_set(expt_in) then begin
  tags=tag_names(expt_in)
;in the current setting it is assumed that the first occurance of _
;seperates the name of the tag in expt and the tag name inside the
;survey  
   for i=0,n_elements(tags)-1 do begin
      tag_sv=strmid(tags(i),0,strpos(tags(i),'_'))
      tag_p=strmid(tags(i),strpos(tags(i),'_')+1)
      if ((not tag_exist(fid_template.expt,tag_sv)) or (not tag_exist(fid_template.expt.(0),tag_p)))  then begin
;         stop
         Error_check=1
         print,'Error: ',tags(i),' is not a valid EXPT input for set_fiducial'
;         tags_valid='sv1_'+tag_names(fid_template.expt.(0))
         print,'Valid inpts are:'
         print,'sv1_'+tag_names(fid_template.expt.(0))
         print,'sv2_'+tag_names(fid_template.expt.(0))
         print,'sv3_'+tag_names(fid_template.expt.(0))
      endif
   endfor
endif


return,error_check
end


