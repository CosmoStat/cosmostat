function tag_check,str,tag,value=value,outvalue=outvalue

;Written by Adam Amara May 2008
;This routine will check a strucutre to see if a specitic tag
;(i.e. structure element). If the tag exists then 1 is returned and
;the value of that tag is stored in an optional output value. This is 
;an extension of the routine tag_exist.pro with the added ability to
;output the tag value.
;
; INPUT PARAMETERS:     
;       str  -  structure variable to search
;       tag  -  tag name to search for, scalar string
;       outvalue - set this keyword to make the function output the
;                  value of the structure
;
; OUTPUTS:
;       Function returns 1b if tag name exists or 0b if it does not.
;
; OPTIONAL OUTPUT KEYWORD:
;       value = the value that the keyword is set to
                             


exist=tag_exist(str,tag,index=index)
if (exist eq 1b) then value=str.(index) else value='NONE'
if keyword_set(outvalue) then exit=value

return,exist
end

