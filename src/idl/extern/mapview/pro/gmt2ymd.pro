
PRO GMT2YMD, GMT, YMD
;+
; NAME:
;   GMT2YMD
; PURPOSE:
;   Convert the YYYYDDD portion of a GMT string to the format YYYY:MM:DD
; EXPLANATION:
;   Converts the YYYYDDD portion of a GMT string to the format
;   YYYY:MM:DD.  The input GMT string can be of any length, but the first 7
;   characters must be in the above form.  The output is a string of length 10.
; CALLING SEQUENCE:
;     GMT2YMD, GMT, YMD
; INPUT:
;     GMT - scalar string containing the first 7 characters in YYYYDDD (Year,
;           Day) format.
; OUTPUT:
;     YMD - 10 character scalar string containing the date in YYYY:MM:DD
;           format.
; EXAMPLE:
;     IDL> GMT2YMD,'20023291849580204716',ymd & print,ymd
;          ===> 2002:11:25
; REVISION HISTORY:
;-
  If N_params() LT 2 then begin
       print,'Syntax - GMT2YMD, GMT, YMD' 
       return
  endif
; Extract the year from the input
 YMD = STRMID(GMT,0,4)
 Iyear = LONG(YMD)
 Data = [1,32,60,91,121,152,182,213,244,274,305,335]
; Test for leap year, ignoring certain century years
 IF( Iyear MOD 4 EQ 0 )THEN Data[2:*] = Data[2:*]+1

; Compute the numerical day and month
 Iday = LONG(STRMID(GMT,4,3))
 Index = MAX(WHERE(Data LE Iday))
 Iday = Iday - Data[Index] + 1
 Imonth = Index+1

; Convert to strings
 Day = STRTRIM(Iday,2)
 IF( STRLEN(Day) EQ 1 )THEN Day = '0' + Day
 Month = STRTRIM(Imonth,2)
 IF( STRLEN(Month) EQ 1 )THEN Month = '0' + Month

; Construct the YMD string
 YMD = YMD + ':' + Month + ':' + Day

 RETURN
 END

