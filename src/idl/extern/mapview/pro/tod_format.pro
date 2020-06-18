FUNCTION TOD_Format, neles, NOHK=nohk, RECSIZE=recsize,SIZE=sz
;+
; NAME:
;	TOD_Format()
; PURPOSE:
;	Defines the time-ordered data (TOD) structure format.
; CALLING SEQUENCE:
;       str = TOD_FORMAT( [Nelements, /NOHK, RECSIZE=,/SIZE )
; INPUTS:
;	neles - The number of elements in the returned array.  Defaults to 1.
; RETURNED:
;	An array of structures.  See the body of this routine for details.
; KEYWORDS:
;	/NOHK     - If present and nonzero, only the time stamp and science
;	           data are defined in the destination structure
;	RECSIZE  - Returns the size of a single archive disk record, in bytes.
;	SIZE     - If present and nonzero, the size of a single archive disk
;	           in bytes is returned instead of a structure.  This overrides
;	           all other keywords.
; NOTES:
;	This function exists to allow define this structure in a single
;	location, allowing easy maintenance of the function.
;
;	The position and velocity data are in celestial coordinates.
;
;	The DISKFMT keyword overrides the NOHK keyword.
;
;	The three types of structures that this routine can return each
;	have different names, so they can coexist.
; MODIFICATION HISTORY:
;	Written by Michael R. Greason, Raytheon STX, 14 January 1998.
;	Nomenclature changes.  DISKFMT keyword added.  MRG, RSTX, 24 April 1998.
;       Nomenclature change---calib structure member becomes gain.  MRG,
;         RSTX, 28 May 1998
;	Pointing changes.  MRG, RSTX, 04 June 1998.
;	Uberframe changes.  MRG, RSTX, 17 September 1998.
;	Real*8 and new flag changes.  MRG, RITSS, 03 July 2001.
;	Corrected the DAFlags data type.  MRG, RITSS, 02 August 2001.
;	Geocentric position and velocity.  MRG, RITSS, 07 August 2001.
;	Calibration replaced by a calibrated data flag.  Science data made
;	  floating point.  MRG, SSAI, 09 May 2002.
;-
on_error, 2
;
If (n_params() LT 1) Then neles = 1
nel = neles > 1
;
recsize = 105842L
If (keyword_set(sz)) Then Return, recsize
;
;			Disk-format version.
;
If (keyword_set(dfmt)) Then Begin
;
;				Instrument Housekeeping - AEU
;
       swp = {AEUSweepDisk,	       $
	       PDU_Analog:intarr(97),  $       ; PDU Analog data block.
	       AEU_Analog1:intarr(57), $       ; AEU Analog data block - 1.
	       AEU_Analog2:intarr(57), $       ; AEU Analog data block - 2.
	       AEU_Window1:intarr(16), $       ; AEU Analog data window - 1.
	       AEU_Window2:intarr(16)}         ; AEU Analog data window - 2.
;
       aeu = {AEUFrameDisk,	       $
	       Day:0, Time:0L,         $       ; Timestamp of AEU frame.
	       Sweep1:swp,	       $       ; First AEU data sweep.
	       Sweep2:swp,	       $       ; Second AEU data sweep.
	       Counters:intarr(5)}	       ; Status counters.
;
;			       Instrument Housekeeping - DEU
;
       deu = {DEUFrameDisk,	       $
	       Day:0, Time:0L,         $       ; Timestamp of DEU frame.
	       Data:intarr(42)         $
	     }
;
;			       Science subframe definition.
;
	If (keyword_set(trnc)) Then Begin
		sci = {SciFrameDiskTr,	  $
			Day:0, Time:0L,	  $     ; Timestamp of sci. frame.
			GenFlags:0L, 	  $     ; Quality flags.
			DAFlags:intarr(10),$	
			Error1:intarr(2), $     ; Telemetry error 1 codes, both packets.
			Error2:intarr(2), $     ; Telemetry error 2 codes, both packets.
			Q:fltarr(4,15)	  $     ; Q-band data.
		      }
	EndIf Else Begin
		sci = {SciFrameDisk,	  $
			Day:0, Time:0L,   $     ; Timestamp of sci. frame.
			GenFlags:0L, 	  $     ; Quality flags.
			DAFlags:intarr(10),$	
			Error1:intarr(2), $     ; Telemetry error 1 codes, both packets.
			Error2:intarr(2), $     ; Telemetry error 2 codes, both packets.
			K:fltarr(4,12),	  $     ; K-band data.
			Ka:fltarr(4,12),  $     ; Ka-band data.
			Q:fltarr(8,15),	  $     ; Q-band data.
			V:fltarr(8,20),	  $     ; V-band data.
			W:fltarr(16,30)	  $     ; W-band data.
		      }
	EndElse
;
;			       The archive record, in on-disk format.
;
	If (keyword_set(trnc)) Then Begin
           ele = {ArchiveDiskTr,       $
	       AEU:aeu, 	       $     ; Instrument housekeeping data.
	       NDEU:0L,  	       $     ; # of valid DEU frames.
	       DEU:replicate(deu, 4),  $
	       Quaternions:fltarr(4,33), $   ; Collection of quaternions.
	       Position:dblarr(3),     $     ; Current WMAP position wrt the Sun.
	       Velocity:dblarr(3),     $     ; Current WMAP velocity wrt the Sun.
	       GeoPosition:dblarr(3),  $     ; Current WMAP position wrt Earth.
	       GeoVelocity:dblarr(3),  $     ; Current WMAP velocity wrt Earth.
	       CalFlag:0,              $     ; Calibrated data flag.
	       NSci:0L,  	       $     ; # of valid science frames.
	       Sci:replicate(sci,30)   $     ; Science data.
	     }
	EndIf Else Begin
           ele = {ArchiveDisk,	       $
	       AEU:aeu, 	       $       ; Instrument housekeeping data.
	       NDEU:0L,  	       $       ; # of valid DEU frames.
	       DEU:replicate(deu, 4),  $
	       Quaternions:fltarr(4,33), $     ; Collection of quaternions.
	       Position:dblarr(3),     $       ; Current MAP position wrt the Sun.
	       Velocity:dblarr(3),     $       ; Current MAP velocity wrt the Sun.
	       GeoPosition:dblarr(3),  $       ; Current MAP position wrt Earth.
	       GeoVelocity:dblarr(3),  $       ; Current MAP velocity wrt Earth.
	       CalFlag:0,              $       ; Calibrated data flag.
	       NSci:0L,  	       $       ; # of valid science frames.
	       Sci:replicate(sci,30)   $       ; Science data.
	     }
	EndElse
;
EndIf Else If (NOT keyword_set(nohk)) Then Begin
;
;				Instrument Housekeeping - AEU
;
       swp = {AEUSweep, 	       $
	       PDU_Analog:lonarr(97),  $       ; PDU Analog data block.
	       AEU_Analog1:lonarr(57), $       ; AEU Analog data block - 1.
	       AEU_Analog2:lonarr(57), $       ; AEU Analog data block - 2.
	       AEU_Window1:lonarr(16), $       ; AEU Analog data window - 1.
	       AEU_Window2:lonarr(16)  $       ; AEU Analog data window - 2.
	     }
       aeu = {AEUFrame, 	       $
	       Day:0, Time:0L,         $       ; Timestamp of AEU frame.
	       Sweep1:swp,	       $       ; First AEU data sweep.
	       Sweep2:swp,	       $       ; Second AEU data sweep.
	       Counters:lonarr(5)      $       ; Status counters.
	     }
;
;			       Instrument Housekeeping - DEU
;
       deu = {DEUFrame, 	       $
	       Day:0, Time:0L,         $       ; Timestamp of DEU frame.
	       Data:lonarr(42)         $
	     }
;
;			       Science subframe definition.
;
       If (keyword_set(trnc)) Then Begin
	       sci = {SciFrameTr, 	 $
		       Day:0, Time:0L,   $     ; Timestamp of sci. frame.
		       GenFlags:0L,	 $     ; Quality flags.
		       DAFlags:lonarr(10),$	
		       Error1:lonarr(2), $     ; Telemetry error 1 codes, both packets.
		       Error2:lonarr(2), $     ; Telemetry error 2 codes, both packets.
		       Q:fltarr(4,15)	 $     ; Q-band data.
		     }
       EndIf Else Begin
	       sci = {SciFrame, 	 $
		       Day:0, Time:0L,   $     ; Timestamp of sci. frame.
		       GenFlags:0L,	 $     ; Quality flags.
		       DAFlags:lonarr(10),$	
		       Error1:lonarr(2), $     ; Telemetry error 1 codes, both packets.
		       Error2:lonarr(2), $     ; Telemetry error 2 codes, both packets.
		       K:fltarr(4,12),   $     ; K-band data.
		       Ka:fltarr(4,12),  $     ; Ka-band data.
		       Q:fltarr(8,15),   $     ; Q-band data.
		       V:fltarr(8,20),   $     ; V-band data.
		       W:fltarr(16,30)   $     ; W-band data.
		     }
       EndElse
;
;			       The archive record, in on-disk format.
;
       If (keyword_set(trnc)) Then Begin
           ele = {ArchiveTr,           $
	       AEU:aeu, 	       $       ; Instrument housekeeping data.
	       NDEU:0L,  	       $       ; # of valid DEU frames.
	       DEU:replicate(deu, 4),  $
	       Quaternions:dblarr(4,33), $    ; Collection of quaternions.
	       Position:dblarr(3),     $      ; Current WMAP position wrt the Sun.
	       Velocity:dblarr(3),     $      ; Current WMAP velocity wrt the Sun.
	       GeoPosition:dblarr(3),  $      ; Current WMAP position wrt Earth.
	       GeoVelocity:dblarr(3),  $      ; Current WMAP velocity wrt Earth.
	       CalFlag:0L,             $      ; Calibrated data flag.
	       NSci:0L,  	       $      ; # of valid science frames.
	       Sci:replicate(sci,30)   $      ; Science data.
	     }
       EndIf Else Begin
           ele = {Archive,             $
	       AEU:aeu, 	       $      ; Instrument housekeeping data.
	       NDEU:0L,  	       $      ; # of valid DEU frames.
	       DEU:replicate(deu, 4),  $
	       Quaternions:dblarr(4,33), $     ; Collection of quaternions.
	       Position:dblarr(3),     $      ; Current WMAP position wrt the Sun.
	       Velocity:dblarr(3),     $      ; Current WMAP velocity wrt the Sun.
	       GeoPosition:dblarr(3),  $      ; Current WMAP position wrt Earth.
	       GeoVelocity:dblarr(3),  $      ; Current WMAP velocity wrt Earth.
	       CalFlag:0L,             $      ; Calibrated data flag.
	       NSci:0L,  	       $      ; # of valid science frames.
	       Sci:replicate(sci,30)   $       ; Science data.
	     }
	EndElse
;
EndIf Else Begin
;
       If (keyword_set(trnc)) Then Begin
	       sci = {SciFrameShortTr, 	 $
		       Day:0, Time:0L,   $     ; Timestamp of sci. frame.
		       GenFlags:0L,	 $     ; Quality flags.
		       DAFlags:lonarr(10),$	
		       Q:fltarr(4,15)	 $     ; Q-band data.
		     }
	       ele = {ArchiveShortTr, NSci:0L, Sci:replicate(sci, 30)}
       EndIf Else Begin
	       sci = {SciFrameShort, 	 $
		       Day:0, Time:0L,   $     ; Timestamp of sci. frame.
		       GenFlags:0L,	 $     ; Quality flags.
		       DAFlags:lonarr(10),$	
		       K:fltarr(4,12),   $     ; K-band data.
		       Ka:fltarr(4,12),  $     ; Ka-band data.
		       Q:fltarr(8,15),   $     ; Q-band data.
		       V:fltarr(8,20),   $     ; V-band data.
		       W:fltarr(16,30)   $     ; W-band data.
		     }
	       ele = {ArchiveShort, NSci:0L, Sci:replicate(sci, 30)}
       EndElse
;
EndElse
;
If (nel GT 1) Then ele = replicate(ele, nel)
;
Return, ele
End
