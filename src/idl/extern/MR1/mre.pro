

!quiet=1
!PATH = EXPAND_PATH("+$CEA_MR_DIR/idl/contrib/astron") + ":" + EXPAND_PATH("+$CEA_MR_DIR/idl/pro") + ":" + !PATH 
 
!help_path= '$CEA_MR_DIR/idl/help:'+!help_path

; resize IDL memory in order to have enough memory for compilation
; .size 65000 20000

; prompt definition

print, " "
print, "*************************************"
print, "** IDL Multiresolution Environment **"
print, "*************************************"
print, " "

print, "HELP: type mrh for the multiresolution help ... "

!PROMPT='MRE> '

astrolib
COMMON XLIVE, XIma

.r alias
