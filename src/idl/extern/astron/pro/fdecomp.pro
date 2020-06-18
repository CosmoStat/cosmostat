pro fdecomp, filename, disk, dir, name, qual, version, OSfamily = osfamily
;+
; NAME:
;     FDECOMP
; PURPOSE:
;     Routine to decompose a file name for any operating system
;
; CALLING SEQUENCE:
;     FDECOMP, filename, disk, dir, name, qual, version, [OSFamily = ]
;
; INPUT:
;     filename - string file name, scalar
;
; OUTPUTS:
;     All the output parameters are scalar strings
;       disk - disk name, always '' on a Unix machine, scalar string
;       dir - directory name, scalar string
;       name - file name, scalar string
;       qual - qualifier, set equal to the characters beyond the last "."
;       version - version number, always '' on a non-VMS machine, scalar string
;
; OPTIONAL INPUT KEYWORD:
;     OSFamily - one of the four scalar strings specifying the operating 
;             system:  'vms','Windows','MacOS' or 'unix'.    If not supplied,
;             then !VERSION.OS_FAMILY is used to determine the OS.
; EXAMPLES:
;     Consider the following file names 
;
;     Unix:    file = '/rsi/idl40/avg.pro' 
;     VMS:     file = '$1$dua5:[rsi.idl40]avg.pro;3
;     Mac:     file = 'Macintosh HD:Programs:avg.pro'
;     Windows: file =  'd:\rsi\idl40\avg.pro'
;       
;     then IDL> FDECOMP,  file, disk, dir, name, qual, version
;       will return the following
;
;                 Disk             Dir          Name        Qual     Version
;       Unix:      ''            '/rsi/idl40/'  'avg'       'pro'       ''
;       VMS:     '$1$dua5'       '[RSI.IDL40]'  'avg'       'pro'       '3'
;       Mac:     'Macintosh HD'  ':Programs:'   'avg'       'pro'       ''
;       Windows:    'd:'         \rsi\idl40\    'avg'       'pro'       ''
;
; NOTES:
;     (1) All tokens are removed between
;           1) name and qual  (i.e period is removed)
;           2) qual and ver   (i.e. VMS semicolon is removed)
;     (2) On VMS the filenames "MOTD" and "MOTD." are distinguished by the 
;         fact that qual = '' for the former and qual = ' ' for the latter.
;
;     A version of FDECOMP that accepts vector input strings is available for
;     IDL V5.3 or later from http://idlastro.gsfc.nasa.gov/ftp/v53/
; ROUTINES CALLED:
;     Function GETTOK()
; HISTORY
;     version 1  D. Lindler  Oct 1986
;     Include VMS DECNET machine name in disk    W. Landsman  HSTX  Feb. 94
;     Converted to Mac IDL, I. Freedman HSTX March 1994          
;     Converted to IDL V5.0   W. Landsman   September 1997
;-
;--------------------------------------------------------
;
  On_error,2                            ;Return to caller

  if N_params() LT 2 then begin
     print, 'Syntax - FDECOMP, filename, disk, [dir, name, qual, ver ] '
     return
  endif

; Find out what machine you're on, and take appropriate action.
 if not keyword_set(OSFAMILY) then osfamily = !VERSION.OS_FAMILY

 case OSFAMILY of

  "MacOS": begin

; disk name is all characters up to the first colon
; directory is string of folders         
; file name+qualifier is all characters after the last colon
; version   is null string
   
  st = filename
  if strpos(st,':') GE 0 then disk = gettok(st,':')  else disk = ''
         
     dir = ':' & tok = ''
     REPEAT BEGIN
        oldtok = tok
        tok = gettok(st,':')
        dir = dir + oldtok + ':'   
     ENDREP UNTIL tok EQ ''

       dir = strmid(dir,1,strpos(dir,oldtok)-1)   
         
     fname = oldtok     & qual = ''
     pos = strpos(fname,'.')
     if pos GE 0 then begin
       name = gettok(fname,'.')
       qual   = fname
     endif 

    version = ''
        
        end

 "vms":  begin                     ; begin VMS version

    st = filename

; get disk

    nodepos = strpos(st,'::')          ; Node name included in directory?
    if nodepos GE 0 then begin
        disk = strmid(st,0,nodepos+2) 
        st = strmid(st,nodepos+2, 999 )
    endif else disk = ''
    if strpos(st,':') GE 0 then disk = disk + gettok(st,':') + ':' else $
                                disk = disk + ''

; get dir

    if strpos( st, ']' ) GE 0 then dir = gettok( st, ']' ) + ']' else dir=''
    if strpos( st, ']' ) GE 0 then dir = dir + gettok( st, ']' ) + ']' 

; get name

    sv_name = st
    name = gettok(st,'.')

; get qualifier

    if (name + '.') EQ sv_name then qual = ' ' else $
    qual = gettok(st,';')

; get version

    version = st

  end   ;  end VMS version

 "Windows": begin

     st = filename
     pos = strpos( st, ':')                 ; DOS diskdrive (i.e. c:)
     if (pos gt 0) then disk = gettok(st,':') + ':' else disk=''

;  Search the path name (i.e. \dos\idl\) and locate all backslashes

     lpos = -1  ; directory position path (i.e. \dos\idl\)
     pos = -1
     repeat begin
        pos = strpos(st, '\',pos+1)
        if (pos GE 0) then lpos = pos
     endrep until pos lt 0

     ;  Parse off the directory path 

     if lpos ge 0 then begin
        dir = strmid(st, 0, lpos+1)
        len = strlen(st)
        if lpos eq (len-1) then $
                st = '' else st = strmid(st,lpos+1,len-lpos-1)
     endif else dir=''

; get Windows name and qualifier (extension)...qual is optional

     lpos=-1
     repeat begin                               
        pos = strpos(st,'.',pos+1)
        if (pos ge 0) then lpos = pos
    endrep until pos lt 0

    ; Parse name and qual (if a qual was found )

     if lpos ge 0 then begin
        len = strlen(st)
        name = strmid(st,0,lpos)
        qual = strmid(st,lpos+1,len-lpos-1)
     endif else begin
        name = st
        qual = '' 
     endelse

     version = ''               ; no version numbers in Windows         
     end

 ELSE: begin

    st = filename

; get disk

    disk = ''

; get dir

    lpos = -1
    pos = -1
    repeat begin
            pos = strpos(st, '/', pos+1)
            if (pos GE 0) then lpos = pos
    endrep until pos LT 0

    if lpos GE 0 then begin
            dir = strmid(st, 0, lpos+1)
            len = strlen(st)
            if lpos eq (len-1) then st = '' else $
                                    st = strmid(st,lpos+1,len-lpos-1)
    endif else dir = ''

; get name and qual

    pos = -1
    lpos = -1
    repeat begin
             pos = strpos(st,'.',pos+1)
             if (pos GE 0) then lpos = pos
    endrep until pos LT 0

    if lpos GE 0 then begin
             len = strlen(st)
             name = strmid(st,0,lpos)
             qual = strmid(st,lpos+1,len-lpos-1)
     endif else begin
         name = st
         qual = '' 
     endelse

    version = ''

 end 

ENDCASE                         ; end OTHER version


  return
  end
