Function DispGMT, gmt, Date=dt, Local=loc
;+
; NAME:
;	DispGMT
; PURPOSE:
;	Converts a time formatted as a MAP GMT into a more
;	human-readable form.
; CALLING SEQUENCE:
;	disp = DispGMT(gmt,[ /DATE, LOCAL= )
; INPUTS:
;	gmt   - The GMT-formatted input time.  This is a scalar string with 
;	        the format:
;	            YYYYDDDHHMMSStttt000
;	        where:	YYYY - the year,
;			DDD  - the day of year,
;			HH   - the time of day, hours,
;			MM   - the time of day, minutes,
;			SS   - the time of day, seconds,
;			tttt - the time of day, tenths of a millisecond.
;			000  - padding out to 20 characters.
; RETURNED:
;	disp  - The string containing the converted date, in GMT epoch:
;		    dd-mon-year hh:mm:ss
; OPTINAL INPUT KEYWORDS:
;	/Date  - If present and nonzero, only the date section is returned.
;	 Local - If present and nonzero, the time is returned as local
;	        time instead of GMT.  The value of local indicates the
;	        timezone to display:
;			 1 = Eastern Standard Time		(EST)
;			-1 = Eastern Daylight Savings Time	(EDT)
;			 2 = Central Standard Time		(CST)
;			-2 = Central Daylight Savings Time	(CDT)
;			 3 = Mountain Standard Time		(MST)
;			-3 = Mountain Daylight Savings Time	(MDT)
;			 4 = Pacific Standard Time		(PST)
;			-4 = Pacific Daylight Savings Time	(PDT)
;	        The timezone abbreviation will be appended to the output string.
; COMMENTS:
;	The GMT is converted into a Julian date, which is then converted
;	into 
; EXAMPLE:
;       IDL> print,dispgmt('20023291849580200000')
;        ===> 25-Nov-2002 18:49:58 GMT
; PROCEDURES USED:
;       JUL2DATE               ;Map Library
;       MONTH_CNV              ;IDL Astro Library
; MODIFICATION HISTORY:
;	Michael R. Greason, Raytheon ITSS, 12 January 1999
;-
on_error, 2
;
;			Check arguments.
;
If (n_params() LT 1) Then message, 'Syntax:  disp = DispGMT(gmt)'
;
;			Constants.
;
If (keyword_set(loc)) Then lc = loc Else lc = 0
Case (lc) Of
	 1 : Begin
	     tzone = ' EST'
	     tzmod = 5
	     End
	-1 : Begin
	     tzone = ' EDT'
	     tzmod = 4
	     End
;
	 2 : Begin
	     tzone = ' CST'
	     tzmod = 6
	     End
	-2 : Begin
	     tzone = ' CDT'
	     tzmod = 5
	     End
;
	 3 : Begin
	     tzone = ' MST'
	     tzmod = 7
	     End
	-3 : Begin
	     tzone = ' MDT'
	     tzmod = 6
	     End
;
	 4 : Begin
	     tzone = ' PST'
	     tzmod = 8
	     End
	-4 : Begin
	     tzone = ' PDT'
	     tzmod = 7
	     End
;
	Else : Begin
	     tzone = ' GMT'
	     tzmod = 0
	     End
EndCase
;
;			Convert the GMT into Gregorian time.
;
jcnv = gmt2jul(gmt) - (double(tzmod) / 24.0d)
Jul2Date, jcnv, year, mon, day, hr
If ((year LT  50) OR ((year LT 200) AND (year GE 100))) Then year = year + 2000
If (year LT 100) Then year = year + 1900
hours   = fix(hr)
minutes = (hr * 60) mod 60
seconds = fix((minutes * 60) mod 60 + 0.5)
minutes = fix(minutes)
;
;			Format the output string.
;
If (keyword_set(dt)) Then Begin
	Disp = string(				$
		day, 				$
		Month_Cnv(mon, /short),	        $
		year, 				$
		Format="(I2.2,'-',A3,'-',I4.4)")
EndIf Else Begin
	Disp = string(				$
		day, 				$
		Month_Cnv(mon, /short),	$
		year, 				$
		hours, 				$
		minutes, 			$
		seconds, 			$
		Format="(I2.2,'-',A3,'-',I4.4,' ',I2.2,':',I2.2,':',I2.2)")
EndElse
;
Return, Disp + tzone
End
