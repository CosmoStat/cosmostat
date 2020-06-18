
 
function format_image, NameImage
Format = 'unknown'
; find the image format
if rstrpos(NameImage, '.d') GT 0  then Format = 'raw' else $
if rstrpos(NameImage, '.fit') GT 0 then Format = 'fits' else $
if rstrpos(NameImage, '.bdf') GT 0 then Format = 'midas' 
return, Format
end 

function rraw, file_name, SIZEX=Sizex, SIZEY=Sizey, TYPE=Type

if keyword_set(SIZEX) then Nc = long(Sizex) else Nc = long(0)  
if keyword_set(SIZEY) then Nl = long(Sizey) else Nl = long(0)
if keyword_set(TYPE) then TypeData = Type else TypeData = ''

if Nl EQ 0 OR NC EQ 0 then BEGIN
   print, "file_name = ", file_name
   p = rstrpos(file_name, ".d")
   if p LT 0 then goto, bad
   f = str_sep(file_name, ".d")
   f = f(0)
   f = str_sep(f, "_")
   s = (size(f))(1)
   if s LT 4 then goto, bad
   Nl = long(f(s-2))
   Nc = long(f(s-1))
   TypeData = f(s-3)
   print, "Nc = ", Nc
   print, "Nl = ", Nl
   print, "type = ", TypeData
   END
taille = Nl*Nc

on_ioerror, bad
openr, unit, file_name, /get_lun

case TypeData of
  'f':  BEGIN 
          a = fltarr (taille)
          b = fltarr (Nc,Nl)
          readu, unit, a
          b(*,*) = a(*)
        END
  'i':  BEGIN
          a = lonarr(taille)
          b = lonarr(Nc, Nl)
          readu, unit, a
          b(*,*) = a(*)
        END
  's':  BEGIN
          a = intarr(taille)
          b = intarr(Nc, Nl)
          readu, unit, a
          b(*,*) = a(*)
        END
  'b':  BEGIN
          a = bytarr(taille)
          b = bytarr(Nc, Nl)
          readu, unit, a
          b(*,*) = a(*)
        END
  'cf': BEGIN
          a = fltarr (2*taille)
          readu, unit, a
          b = complexarr(Nc, Nl)
          i = long(0)
          j = long(0)
          N = long(Nc)
          for i = 0, Nl-1 do begin
          for j = 0, Nc-1 do begin
             b(i,j) = complex(a(2*i*N+j), a(2*i*N+j+1))
          endfor
          endfor

        END
  'd':  BEGIN
          a = dblarr (taille)
          readu, unit, a
          b = dblarr(Nc, Nl)
          b(*,*) = a(*)
        END 
  'cd': BEGIN
          a = dblarr (2*taille)
          b = complexarr(Nc, Nl)
          readu, unit, a
          i = long(0)
          j = long(0)
          N = long(Nc)
          for i = 0, Nl-1 do begin
          for j = 0, Nc-1 do begin
             b(i,j) = complex(a(2*i*N+j), a(2*i*N+j+1))
          endfor
          endfor
        END 
   ELSE: goto, bad
endcase

free_lun, unit
goto, done
bad: print, 'Error: Unable to open file: ', file_name, ' ',!err_string
b = -1
done: 
return, b
end
;--------------------------------------------------------------

function rim, NameImage
if N_PARAMS() LT 1 then BEGIN
               print, 'result = rim(NameImage)' 
               return, -1
               END
Format = format_image(NameImage)
CASE Format OF
   'midas': BEGIN
            File = str_sep(NameImage,'.bdf')
            File = File(0)
            MID_RD_IMAGE, File, Dat, NAXIS, NPIX
            END
   'fits': Dat = readfits(NameImage)
   'raw': Dat = rraw(NameImage)
    else: BEGIN
          print, 'rim: unable to read data, unknown format'
          Dat = -1
          END
ENDCASE
return, Dat
END
;--------------------------------------------------------------

function def_name_raw, NameStep, Dat
Nl = (size(Dat))(2)
Nc = (size(Dat))(1)
p = rstrpos(NameStep, ".d")
file = NameStep
type = 'f'
if p LT 0 then $
BEGIN
  case (size(Dat))(3) of
       0:  goto, bad
       1:  type='b'
       2:  type='s'
       3:  type='i'
       4:  type='f'
       5:  type='d'
       6:  type='cf'
      ELSE: goto, bad
   endcase
  
   file = file + '_' + type + '_' +  string(Nl) + '_' + string(Nc) + '.d'
   file = strcompress(file, /REMOVE_ALL)
END
print, file
goto, DONE

BAD: file=-1
     
DONE: 
return, file
END

;------------------------------------------------------

pro wraw, file_name, b, OUTPUTFILE=OutputFile
print, "file name = ", file_name
Nl = (size(b))(2)
Nc = (size(b))(1)
print, "Nl = ", Nl
print, "Nc = ", Nc
taille= long(nl)*long(nc)
print, "Size = ", taille
p = rstrpos(file_name, ".d")
file = file_name
type = 'f'
case (size(b))(3) of
       0: goto, bad
       1: BEGIN
            if p LT 0 then type='b'
            a = bytarr(taille)
            a(*) = b(*,*)
          END
       2: BEGIN
            if p LT 0 then type='s'
            a = intarr(taille)
            a(*) = b(*,*)
          END
       3: BEGIN
            if p LT 0 then type='i'
            a = lonarr(taille)
            a(*) = b(*,*)
          END
       4: BEGIN
           if p LT 0 then type='f'
           a = fltarr(taille)
           a(*) = b(*,*)
          END
       5: BEGIN
           if p LT 0 then type='d'
           a = dblarr(taille)
           a(*) = b(*,*)
          END
       6: BEGIN
           if p LT 0 then type='cf'
           a = dblarr(2*taille)
           i = long(0)
           j = long(0)
           N = long(Nc)
           for i = 0, Nl-1 do begin
           for j = 0, Nc-1 do begin
             a(2*i+j) = float(b(i,j))
             a(2*i+j+1) = imaginary(b(i,j))
           endfor
          endfor
          END
    ELSE: goto, bad
endcase
if p LT 0 then BEGIN
          file = file + '_' + type + '_' +  string(Nl) + '_' + string(Nc) + '.d'
          file = strcompress(file, /REMOVE_ALL)
          END
print, file
if keyword_set(OUTPUTFILE) then OutputFile = file
openw, unit, file, /get_lun
on_ioerror, bad
writeu, unit, a

goto, done
bad: print, !err_string
done: free_lun, unit

return
end

;-----------------------------------------------------------------

pro wim, NameImage, Dat, FORMAT=Format
if N_PARAMS() LT 2 then BEGIN
               print, 'wim, NameImage, Dat, FORMAT=Format' 
               print, 'format = ''midas'', ''fits'', or ''raw'''
               return 
               END

if NOT keyword_set(FORMAT) then Format = format_image(NameImage)
CASE Format OF
   'midas': BEGIN
            NAXIS=2
            NPIX=lonarr(2)
            NPIX(0) = (size(Dat))(1)
            NPIX(1) = (size(Dat))(2)
            File = str_sep(NameImage,'.bdf')
            File = File(0)
            MID_UP_IMAGE, File, Dat, NAXIS, NPIX
            END
   'fits': BEGIN
           File = NameImage
           p = rstrpos(File, ".fit")
           if p LT 0 then File = File + '.fits'
           writefits, File, Dat
           END
   'raw': wraw, NameImage, Dat
    else: print, 'Error: unknown format data
ENDCASE
END

