; $Id: interpolate_quaternions.pro,v 1.5 2003/11/12 19:50:07 bhill Exp $
;
Pro Interpolate_Quaternions, input_q, offset, interp_q, status, $
	ExtLimit=elim, Transpose=tflg, Dmax=delta_max, Verbose=vflg
;+
; NAME:
;	Interpolate_Quaternions
; PURPOSE:
;	Routine to interpolate quaternions.
; CALLING SEQUENCE:
;	Interpolate_Quaternions, input_q, offset, interp_q [, status]
; INPUTS:
;	input_q  - Set of 4 evenly-spaced quaternions (in a 4x4 array).
;		   See the COMMENTS section for how this array should
;		   be arranged.
;	offset   - Dimensionless time offset relative to the first quaternion.
; OUTPUTS:
;	interp_q - The interpolated quaternion.
;	status   - A status code:
;			 0=normal interpolation,
;			-1=normal extrapolation--low,
;			-2=normal extrapolation--high,
;			 1=offset out of bounds--low,
;			 2=offset out of bounds--high,
; KEYWORDS:
;	ExtLimit  - The dimensionless extrapolation limit; it defines the
;		    permitted range in Offset.  If this optional parameter
;		    is not supplied then it defaults to 1.0, resulting in
;		    a permitted range in Offset of -1 to 4 (0 to 3 indicates
;		    interpolation while anything outside that range indicates
;		    extrapolation).
;	Transpose - If present and nonzero, the transpose of the input_q
;		    array is used.
;	Dmax      - Sets the maximum allowed displacement between adjacent q's.
;                   Physically, |delta_q|^2 = |omega*delta_t/2|^2, the default
;                   is 0.02 rad^2 (omega ~ 10 deg/sec, delta_t=1.536 sec)
;                   which should cover all map applications.  Any displacements
;                   larger than Dmax will produce an error condition.
;       Verbose   - If present and nonzero, the input and intermediate 
;                   quaternions and the result will be printed.
; COMMENTS:
;	This routine expects a unifomly sampled set of quaternions Q1,Q2,Q3,Q4.
;	It interpolate a quaternion for any time between Q1 and Q4, inclusive.
;	The output is calculated at a time T_Out, expressed in terms of the
;	sampling of the input quaternions:
;
;			   T_Out - T(Q1)
;		Offset = -----------------
;			   T(Q2) - T(Q1)
;
;	where T(Q1) is the time at quaternion Q1, and so forth.  That is,
;	the time for the output quaternion (variable OFFSET) should be
;	a number in the range -1.000 to 4.000 inclusive.  Input values outside
;	that range result in an error.  Input values outside 0.0 to 3.0 result
;	in extrapolation instead of interpolation.
;
;       In other words, Offset is essentially a floating point subscript,
;       similar to the those used by the IDL intrinsic routine INTERPOLATE.
;
;	For optimal results, OFFSET should be in the range [1.0, 2.0] -- that
;	is, the input quaternions Q1...Q4 should be arranged such that 2 come
;	before the desired output and 2 come after.
;
;	This routine expects input_q to be ordered as follows, where
;       Q1-Q4 are the 4 input quaternions:
;		input_q = [[Q1], [Q2], [Q3], [Q4]]
;       i.e., the Y subscript selects a quaternion and the X subscript
;       selects an element.
;
;	If passing
;		input_q = transpose([[Q1], [Q2], [Q3], [Q4]])
;	i.e., so that the X subscript selects a quaternion, and the
;       Y subscript selects an element, then use the Transpose keyword.
;
;	All computations are performed in double precision, and a double
;	precision quaternion is returned regardless of the data type of
;	the input quaternions.
; MODIFICATION HISTORY:
;	Adapted from the Fortran by Michael R. Greason, Raytheon ITSS,
;	    03 August 1999.
;	Check of rate of change added in order to interpolate over
;	    apparent discontinuities resulting from use of principal
;	    value only.  G. Hinshaw, July, 2001.
;	Row/column usage made consistent.  RSH, RITSS, 3 Aug 2001.
;	Offset check corrected.  MRG, SSAI, 13 February 2002.
;       Row/column doc corrected.  RSH, 12 Nov 2003.
;-
on_error, 2
;
;		Check arguments.
;
;			Present?
;
If (n_params() LT 3) Then message, $
	'Syntax: Interpolate_Quaternions, input_q, offset, interp_q [, status]'
;
If (n_elements(elim) LE 0) Then extlim = 1.0d0 $
                           Else extlim = double(elim)
;
;			Verify that the offset is within limits.
;
Case (1) Of
	(offset LT (-extlim))		: status =  1
	(offset GT (3.0d0 + extlim))	: status =  2
	(offset LT 0.0d0)		: status = -1
	(offset GT 3.0d0)		: status = -2
	Else				: status =  0
EndCase
If (status GT 0) Then Return

;
;		Set the maximum allowed displacement (0.5*omega*delta_t)^2
;		Default: omega = 10 deg/sec, delta_t = 1.536 sec
;
If ( not keyword_set(delta_max) ) Then delta_max = 0.02
;
;		Compute the weights for interpolation.
;
xp0 = offset - 1.0d0
xn0 = -xp0
xp1 = xp0 + 1.0d0
xn1 = xp0 - 1.0d0
xn2 = xp0 - 2.0d0
w = [	(xn0 * xn1 * xn2 / 6.0d0),	$
	(xp1 * xn1 * xn2 / 2.0d0),	$
	(xp1 * xn0 * xn2 / 2.0d0),	$
	(xp1 * xp0 * xn1 / 6.0d0)	]
;w = [[w], [w], [w], [w]]
;
;		Transpose input, if necessary, and check for sign changes
;
If (keyword_set(tflg)) Then iq = transpose(double(input_q)) $
                       Else iq = double(input_q)
dq2 = total((iq - shift(iq,0,1))^2,1)
For i = 1,3 Do If (dq2[i] gt delta_max) Then iq[*,i:*] = -iq[*,i:*]
;
;		Test again for continuity
;
dq2 = total((iq - shift(iq,0,1))^2,1)
If (max(dq2[1:*]) gt delta_max) Then Begin
    Status=-3
    Return
Endif
;
;		Interpolate within current quaternion array.
;
interp_q = dblarr(4)
For i = 0, 3 Do Begin
    qt = reform(iq[i,*])
    interp_q[i] = total(qt * w)
EndFor
;
;		Normalize the interpolated quaternion.
;
norm = sqrt(total(interp_q ^ 2))
If (norm NE 0) Then interp_q = interp_q / norm
If keyword_set(vflg) Then Begin
    Print,' '
    Print,' Input: '
    Print,  Input_q
    Print,' Inter: '
    Print,  iq
    Print,' Output:', interp_q
Endif

;
Return
End
