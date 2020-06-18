pro uvbybeta,by,m1,c1,beta,n,name,Te,MV,eby,delm0,radius,TEXTOUT=textout
;+
; NAME:
;	UVBYBETA
; PURPOSE:
;	Derive dereddened colors, metallicity, and Teff from Stromgren colors.
; EXPLANATION:
;	Adapted from FORTRAN routine of same name
;	published by T.T. Moon, Communications of University of London
;	Observatory, No. 78.  Can be used either interactively or called
;	from a main procedure.
;
; CALLING SEQUENCE:
;	uvbybeta                    ;Prompt for all parameters
;	uvbybeta,by,m1,c1,beta,n    ;Supply inputs, print outputs
;	uvbybeta, by, m1, c1, beta, n, name, Te, Mv, Eby, delm0, radius, 
;			[ TEXTOUT= ]
;
; INPUTS:
;	by - Stromgren b-y color, scalar
;	m1 - Stromgren line-blanketing parameter, scalar
;	c1 - Stromgren Balmer discontinuity parameter, scalar
;	beta - H-beta line strength index.  If beta is not know UVBYBETA
;		will compute a value based on by, m1,and c1.
;	n -  Integer (1-8) giving approximate stellar classification
;
;	(1) B0 - A0, classes III - V, 2.59 < BETA < 2.88,-0.20 <   c0  < 1.00
;	(2) B0 - A0, class   Ia     , 2.52 < BETA < 2.59,-0.15 <   c0  < 0.40
;	(3) B0 - A0, class   Ib     , 2.56 < BETA < 2.61,-0.10 <   c0  < 0.50
;	(4) B0 - A0, class   II     , 2.58 < BETA < 2.63,-0.10 <   c0  < 0.10
;	(5) A0 - A3, classes III - V, 2.87 < BETA < 2.93,-0.01 < (b-y)o< 0.06
;	(6) A3 - F0, classes III - V, 2.72 < BETA < 2.88, 0.05 < (b-y)o< 0.22
;	(7) F1 - G2, classes III - V, 2.60 < BETA < 2.72, 0.22 < (b-y)o< 0.39
;	(8) G2 - M2, classes  IV _ V, 0.20 < m0   < 0.76, 0.39 < (b-y)o< 1.00
;
;	name - scalar string giving name of star.  Used only when writing to 
;		disk for identification purposes.
;
; OPTIONAL INPUT KEYWORD:
;	TEXTOUT   Used to determine output device.  If not present, the
;	value of !TEXTOUT system variable is used (see TEXTOPEN)
;		textout=1	Terminal with /MORE
;		textout=2	Terminal without /MORE
;		textout=3	uvbybeta.prt   (output file)
;		textout=4	Laser Printer 
;		textout=5      User must open file         
;		textout=7	Append to existing uvbybeta.prt file
;		textout = filename (default extension of .prt)
;
; OPTIONAL OUTPUTS:
;	Te - approximate effective temperature
;	MV - absolute visible magnitude
;	Eby - Color excess b-y
;	delm0 - metallicity index, delta m0, may not be calculable for early
;		B stars.
;	radius - Stellar radius (R/R(solar))
;
; SYSTEM VARIABLES:
;	If keyword textout not used, the non-standard system variable !TEXTOUT 
;	becomes the output device indicator.
;	Set  !TEXTOUT =3 to have results directed to a file UVBYBETA.PRT 
;	If all output parameters were supplied, then type TEXTCLOSE to close
;	this file
;
; REVISION HISTORY:                                           
;	W. Landsman          IDL coding              February, 1988
;	Keyword textout added, J. Isensee, July, 1990
;	Made some constants floating point.   W. Landsman    April, 1994
;	Converted to IDL V5.0   W. Landsman   September 1997
;-
 npar = N_params()

 if not keyword_set( TEXTOUT ) then textout = !TEXTOUT  ;default output dev.

 Rm1 = -0.33 & Rc1 = 0.19 & Rub = 1.53  	;Parameter values
 init = 0

 READ_PAR:  if ( npar LT 4 ) then begin 
  ans = ''
  print,'Enter (b-y), m1, c1, and Beta in that order ([RETURN] to exit)'
  read,'', ans
  if ans eq '' then begin               ;Normal Exit
    if ( init EQ 1 ) then textclose, TEXTOUT = textout
    return 
  endif else ans = getopt(ans)
  if ( N_elements(ans) NE 4 ) then begin
    message, 'INPUT ERROR - Expecting 4 scalar values', /CON
    print, 'Enter 0.0 for beta if it is not known: '
    goto, READ_PAR 
  endif else begin
    by = ans[0] & m1 = ans[1] & c1 = ans[2]  & beta = ans[3]
    endelse
 endif

 flag1 = (beta EQ 0.)
 flag2 = 0
 ub = c1 + 2*(m1+by)

 READ_GROUP:  if ( npar LT 5 )then begin

   print,' The following group of stars are available'
   print, $ 
     '(1) B0 - A0, classes III - V, 2.59 < BETA < 2.88,-0.20 <   c0  < 1.00'
   print, $
     '(2) B0 - A0, class   Ia     , 2.52 < BETA < 2.59,-0.15 <   c0  < 0.40'
   print, $
     '(3) B0 - A0, class   Ib     , 2.56 < BETA < 2.61,-0.10 <   c0  < 0.50'
   print, $ 
     '(4) B0 - A0, class   II     , 2.58 < BETA < 2.63,-0.10 <   c0  < 0.10'
   print, $ 
     '(5) A0 - A3, classes III - V, 2.87 < BETA < 2.93,-0.01 < (b-y)o< 0.06'
   print, $
     '(6) A3 - F0, classes III - V, 2.72 < BETA < 2.88, 0.05 < (b-y)o< 0.22'
   print,$ 
     '(7) F1 - G2, classes III - V, 2.60 < BETA < 2.72, 0.22 < (b-y)o< 0.39'
   print, $ 
     '(8) G2 - M2, classes  IV _ V, 0.20 < m0   < 0.76, 0.39 < (b-y)o< 1.00
   n = 0                   
   read,'Enter group number to which star belongs: ',n

 endif

 case n of

 1: BEGIN
    Eby = ( 13.608*by-ub+1.467 ) / (13.608-Rub)
    DEREDD, Eby, by, m1, c1, ub, by0, m0, c0, ub0
    IF flag1 EQ 1 then beta = -0.027352*c0^3 + 0.161463*c0^2 $
                             + 0.132557*c0   + 2.61033
   g = ALOG10(beta - 2.515) - 1.6*ALOG10(c0 +0.322)
   MV = 3.4994 + 7.2026*ALOG10(beta - 2.515) -2.3192*g + 2.9375*g^3
   Te = 5040/(0.2917*c0 + 0.2)
   m0zams = 0.0957758*c0^3-0.139003*c0^2+0.109804*c0 + 0.07473
   delm0 = m0zams - m0
   flag2 = 1
   END

 2: BEGIN
    Eub = ( 1.5*c1 - ub + 0.035) / (1.5/(Rub/Rc1)-1)
    Eby = Eub/Rub
    DEREDD, Eby, by, m1, c1, ub, by0, m0, c0, ub0
    if ( flag1 EQ 1 ) then beta = 0.037*c0 + 2.542
    END

 3: BEGIN
    Eub = (1.36*c1-ub+0.004) / (1.36/(Rub/Rc1)-1)
    Eby = Eub/Rub
    DEREDD, Eby, by, m1, c1, ub, by0, m0, c0, ub0
    if flag1 eq 1 then beta = 0.047*c0 +2.578
    END

 4: BEGIN
    Eub = ( 1.32*c1 - ub - 0.056) / ( 1.32 / (Rub/Rc1)-1 )
    Eby = Eub/Rub
    DEREDD, Eby, by, m1, c1, ub, by0, m0, c0, ub0
    if ( flag1 EQ 1 ) then beta = 0.066*c0+2.59
    END

 5: BEGIN
    m = m1 - Rm1*by
    by0 = 4.2608*m^2 - 0.53921*m - 0.0235
    REPEAT BEGIN
       bycorr = by0
       m0 = m1 - Rm1*(by-bycorr)
       by0 = 14.0881*m0^2 - 3.36225*m0 + 0.175709
    ENDREP UNTIL ( abs(bycorr - by0) LT 0.001)
    Eby = by - by0
    DEREDD, Eby, by, m1, c1, ub, by0, m0, c0, ub0
    if flag1 eq 1 then beta = 2.7905 - 0.6105*by + 0.5*m0 + 0.0355*c0
    r = 0.35*(c1-Rc1*by) - (beta-2.565)
    a0 = by0+ 0.18*(ub0-1.36)
    MV = 1.5 + 6.0*a0 - 17.0*r
    Te = 5040/( 0.8572*a0 + 0.5152 )
    m0zams = -3.95105*by0^2 + 0.86888*by0 + 0.1598
    delm0 = m0zams - m0
   end

 6: begin
    if flag1 then begin
	print,' Estimate of Beta only valid if star is unreddened'
        beta = 3.06 - 1.221*by - 0.104*c1
    endif
    m1zams = -2.158*beta^2 +12.26*beta-17.209
    if ( beta LE 2.74 ) then begin

	c1zams = 3.0*beta - 7.56
	MVzams = 22.14 - 7*beta

   endif else if ( ( beta GT 2.74 ) and ( beta LE 2.82 ) ) then begin

	c1zams = 2.0*beta - 4.82
	MVzams = 11.16-3*beta

   endif else begin
	c1zams = 2.0*beta-4.83
	MVzams =-88.4*beta^2+497.2*beta-696.41

   endelse        
   delm1 = m1zams - m1
   delc1 = c1-c1zams
   if delm1 lt 0. then $
	by0 = 2.946 - beta - 0.1*delc1 - 0.25*delm1 else $
        by0 = 2.946 - beta - 0.1*delc1
   Eby = by - by0
   Deredd, eby, by, m1, c1, ub, by0, m0, c0, ub0
   delm0 = m1zams - m0
   delc0 = c0 - c1zams
   MV = MVzams -9.0*delc0
   Te = 5040 / (0.771453*by0 + 0.546544)
 end

 7: begin
   if flag1 then begin 
	byinit = by
	m1init = m1
	for i = 1,10 do begin
          m1by = 2.5*byinit^2 - 1.32*byinit + 0.345
          bycorr = byinit + (m1by-m1init) / 2.0
          if ( abs(bycorr-byinit) LE 0.001 ) then goto,T71
	  byinit = bycorr
	  m1init = m1by
	endfor
        T71: beta = 1.01425*bycorr^2 - 1.32861*bycorr + 2.96618  
    endif

    m1zams = 6.41701*beta^2 - 34.4538*beta + 46.4167
    MVzams = 5.48012*beta^3 + 11.0494*beta^2 - 188.748*beta + 324.482

    if beta le 2.65 then $
	c1zams = 2*beta - 4.91 else $
	c1zams = 11.1555*beta^2-56.9164*beta+72.879

     delm1 = m1zams - m1
     delc1 = c1 - c1zams
     dbeta = 2.72 - beta
     by0 = 0.222+1.11*dbeta +2.7*dbeta^2-0.05*delc1-(0.1+3.6*dbeta)*delm1
     Eby = by - by0
     Deredd,Eby,by,m1,c1,ub,by0,m0,c0,ub0
     delm0 = m1zams - m0
     delc0 = c0 - c1zams
     f = 9.0 + 20.0*dbeta
     MV = MVzams - f*delc0
     Te = 5040/(0.771453*by0 + 0.546544)
 end

 8:   begin
     if ( flag1 EQ 1 ) then flag1 = 2
     if ( by LE 0.65 ) then $
           Eby = (5.8651*by - ub -0.8975) / (5.8651 - Rub) $

     else if ( ( by GT 0.65 ) and ( by LT 0.79 ) ) then begin 

           Eby = (-0.7875*by - c1 +0.6585)/(-0.7875 - Rc1)
           by0 = by -Eby
        if ( by0 LT 0.65 ) then $
           Eby = (5.8651*by - ub -0.8975) / (5.8651-Rub)

     endif else begin 

        Eby = ( 0.5126*by - c1 - 0.3645 ) / (0.5126-Rc1)
        by0 = by - Eby
        if ( by0 LT 0.79 ) then $ 
                  Eby = (-0.7875*by - c1 + 0.6585) / (-0.7875-Rc1)
        by0  = by - Eby
        if ( by0 LT 0.65 ) then $ 
                  Eby = ( 5.8651*by - ub - 0.8975) / (5.8651-Rub)

     endelse 

        DEREDD,Eby,by,m1,c1,ub,by0,m0,c0,ub0
        m1zams = 42.93678*by0^4 - 122.466*by0^3 + 122.1875*by0^2 $
               - 49.43695*by0 + 7.18436
        IF by0 lt 0.65 THEN BEGIN
		c1zams = -28.7056*by0^3 +42.7486*by0^2 -21.278*by0 + 3.78514
                MVzams = -552.48*by0^4 + 1272.503*by0^3-1101.257*by0^2 $
                       +  432.156*by0 - 59.2095
        ENDIF ELSE IF (by0 GE 0.65) and (by0 lt 0.79) THEN BEGIN
		c1zams = -0.631821*by0^2+0.116031*by0+0.33657
                MVzams = 1.37632*by0^2 + 4.97911*by0+3.4305
        ENDIF ELSE BEGIN
                c1zams = -0.010028*by0^2 + 0.530426*by0 - 0.37237
                MVzams =  1.18298*by0^2  + 3.92776*by0 + 4.37507
        ENDELSE
	delm0 = m1zams - m0
	delc0 =c0 - c1zams
	IF (by0 LE 0.505) THEN BEGIN
		f = 10. - 80.*(by0-0.38)
                Te = 10^(-0.416*by0+3.924)
        ENDIF ELSE BEGIN
                f = 0.0
                Te = 10^(-0.341*by0+3.869)
        ENDELSE
        MV = MVzams - f*delc0 + 3.2*delm0 - 0.07
      END 
 ELSE: BEGIN
      print,'A stellar group of',n,' is not available
      npar = npar<4
      goto, READ_GROUP 
      end

 endcase
 if (n GE 2) and ( n LE 4 ) then begin
     betaza = -0.160402*c0^4 + 0.277363*c0^3 - 0.099623*c0^2 + $
               0.228638*c0 + 2.62745
     B = betaza - 2.5
     MVzams =203.704*B^3 - 206.98*B^2 + 77.18*b - 9.563
     dbeta = betaza - beta
     dMV = -121.6*dbeta^2 +61.0*dbeta + 0.08
     MV = MVzams - dMV
     Te = 5040 / (0.35866*ub0 + 0.27346)
     flag2 = 2
endif
 if ( by0 LE 0.335 ) then $
            FV = -6.759*by0^3 + 3.731*by0^2 - 1.092*by0 + 3.981 $
       else FV = -0.534*by0 + 3.959
 radius = 10^(2.*(4.236-0.1*MV - FV))
 if ( npar GE 7 ) then return
 if ( flag2 EQ 2 )then metal = 'no delta(m0)' else metal = 'delta(m0) = '
 beta = fix(beta*1000.+0.5)/1000.
 Teff = fix( fix(Te/10.+0.5)*10.)
 if !TEXTUNIT eq 0 then textopen,'uvbybeta',textout=textout
 init = 1                          ;First star has been done
 if (TEXTOUT ne 1) and (npar lt 6) then begin ;Prompt for star name?
     name = ''
     read,'Enter name of star: ',name
 endif
 if (textout NE 1) or (npar GE 6 ) then printf,!TEXTUNIT, $
   form="(10X,'Star is : ',A20,10X,'Processed in group',I2,/)",name,n
 fmt = "(1x,' b-y   = ',F6.3,7X,'m1 = ',F5.3,10X,'c1 = ',F5.3"
 case flag1 of 
    2: printf, !TEXTUNIT, F = fmt + ",6X,'Beta is not used')", by, m1, c1
    1: printf, !TEXTUNIT, $
       f = FMT + ",2x,'estimated Beta  = ',F5.3)", by, m1, c1,beta
    0: printf,!TEXTUNIT, f = FMT + ",8x,'Beta  = ',F5.3)",by,m1,c1,beta
 endcase

 fmt = "(1x,'(b-y)0 = ',F6.3,7X,'m0 = ',F5.3,10X,'c0 = ',F5.3,"
 printf,!TEXTUNIT,form=fmt + "8X,'E(b-y) = ',F6.3,/)",by0,m0,c0,Eby

 printf,!TEXTUNIT,form="(1X,'Absolute Magnitude (Mv) = ',F6.2,5x," + $
       "'Radius  (R/R[solar]) = ',F7.2)",MV,radius

 fmt1 = "(1X,A12,25X,'Effective Temperature (Teff) = ',I5,1X,'K'//)"
 fmt2 = "(1X,A12,F6.3,20X,'Effective Temperature (Teff) = ',I5,1X,'K'//)"

 if ( flag2 EQ 2 ) then printf,!TEXTUNIT,form=fmt1,metal,Teff else  $
                       printf,!TEXTUNIT,form=fmt2,metal,delm0,Teff

 if ( npar LT 4 ) then goto, READ_PAR 
 return
 end
