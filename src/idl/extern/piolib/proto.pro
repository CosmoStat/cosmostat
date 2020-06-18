FUNCTION PIOGETVERSION,VersionTag
;+
; VersionTag : =>  Input: name of the piolib version tag 

;-
FUNCTION PIOCREATEREPOSITORY,DataRepository,HTMLRepository
;+
; DataRepository : =>  Input: name of the data database 
; HTMLRepository : =>  Input: name of the HTML database 

;-
FUNCTION PIOGETBEGINOBJECTIDX,ObjectName,MyGroup
;+
; ObjectName : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOGETENDOBJECTIDX,ObjectName,MyGroup
;+
; ObjectName : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOTYPESIZE,typeIn
;+
; typeIn : =>  Input: the piolib type 

;-
FUNCTION PIOGETVALIDX,Value,Obj,MyGroup
;+
; Value : =>  Input: Value 
; Obj : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOCREATETOIGRP,Groupname,ROIname,Init_Size
;+
; Groupname : =>  Input: name of the opened group 
; ROIname : =>  Input: Ring Information group name 
; Init_Size : =>  Input: default size 

;-
FUNCTION PIOCREATETOIOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENTOIGRP,GroupName,mode
;+
; GroupName : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITETOIOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETOIOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETOIOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETOIOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSETOIGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOOPENTOIOBJECT,MyObjectName,typein,Command,step,mode,MyGroup
;+
; MyObjectName : =>  Input: Object name 
; typein : =>  Input: output data type 
; Command : =>  Input: output data type 
; step : =>  Input: step of reading 
; mode : =>  Input: access mode ("r" | "w" ) 
; MyGroup : =>  Input: link with the opened group 

;-
FUNCTION PIOGETTOIOBJECT,data,MyObject
;+
; data : =>  Input: return table 
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOFLUSHTOIOBJECT,MyObject
;+
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOSETTOIOBJECT,TabIn,MyObject
;+
; TabIn : =>  Input: Input table 
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOSEEKTOIOBJECT,SampleNumber,MyObject
;+
; SampleNumber : =>  Input: Number of sample to skip 
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOCLOSETOIOBJECT,MyObject
;+
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOWRITERINGINDEXGRP,BeginIndx,EndIndx,RingNum,NbNum,MyGroup
;+
; BeginIndx : =>  Input: Begin index table 
; EndIndx : =>  Input: End index table 
; RingNum : =>  Input: ring rank index 
; NbNum : =>  Input: number of ring index 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOGETTOIRINGINDEXGRP,Index,MyGroup
;+
; Index : =>  Input: Absolut index value 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOREADRINGINDEXGRP,BeginIndx,EndIndx,RingNum,MyGroup
;+
; BeginIndx : =>  Output: Begin index table 
; EndIndx : =>  Output: End index table 
; RingNum : =>  Output: ring rank index 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOINFOTOIGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,TOItype=TOItype,TOIname=TOIname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbTOI=NbTOI,ROIGroup=ROIGroup
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; TOItype : =>  Output: TOI format 
; TOIname : =>  Output: TOI name list 
; BeginIndx : =>  Output: first index 
; EndIndx : =>  Output: last index 
; NbTOI : =>  Output: Number of TOI object 
; ROIGroup : =>  Output: ROI group name
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOFREEINFOTOIGRP,FLGname,TOItype,TOIname,BeginIndx,EndIndx
;+
; FLGname : =>  Input: FLAG name list 
; TOItype : =>  Input: TOI format 
; TOIname : =>  Input: TOI name list 
; BeginIndx : =>  Input: first index 
; EndIndx : =>  Input: last index 

;-
FUNCTION PIOCREATETIMESTAMPGRP,Groupname,list_TRII
;+
; Groupname : =>  Input: name of group
; list_TRII : =>  Input: list of TRIIs of group

;-
FUNCTION PIOCREATETIMESTAMPOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENTIMESTAMPGRP,GroupName,mode
;+
; GroupName : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITETIMESTAMPOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETIMESTAMPOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETIMESTAMPOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETIMESTAMPOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSETIMESTAMPGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOOPENTIMESTAMPOBJECT,MyObjectName,type,Command,step,mode,MyGroup
;+
; MyObjectName : =>  Input: Object name 
; type : =>  Input: output data type 
; Command : =>  Input: output data type 
; step : =>  Input: step of reading 
; mode : =>  Input: access mode ("r" | "w" ) 
; MyGroup : =>  Input: link with the opened group 

;-
FUNCTION PIOGETTIMESTAMPOBJECT,data,MyObject
;+
; data : =>  Input: return table 
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOFLUSHTIMESTAMPOBJECT,MyObject
;+
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOSETTIMESTAMPOBJECT,TabIn,MyObject
;+
; TabIn : =>  Input: Input table 
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOSEEKTIMESTAMPOBJECT,SampleNumber,MyObject
;+
; SampleNumber : =>  Input: Number of sample to skip 
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOCLOSETIMESTAMPOBJECT,MyObject
;+
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOWRITEPACKETINDEX,data_timesec,data_period,data_ndata,Ring_num,NbNum,MyGroup
;+
; data_timesec : =>  Input: input TIMESEC table 
; data_period : =>  Input: input PERIOD table 
; data_ndata : =>  Input: input NDATA table 
; Ring_num : =>  Input: list of information index 
; NbNum : =>  Input: number of information 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOGETTIMESTAMPPACKETINDEXGRP,Index,MyGroup
;+
; Index : =>  Input: input index value 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOREADPACKETINDEX,data_timesec,data_period,data_ndata,MyGroup
;+
; data_timesec : =>  Output: Begin index table 
; data_period : =>  Output: End index table 
; data_ndata : =>  Output: End index table 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOINFOTIMESTAMPGRP,MyGroup,TIMESTAMPtype=TIMESTAMPtype,TIMESTAMPname=TIMESTAMPname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbTIMESTAMP=NbTIMESTAMP,PacketGroup=PacketGroup
;+
; TIMESTAMPtype : =>  Output: TIMESTAMP format 
; TIMESTAMPname : =>  Output: TIMESTAMP name list 
; BeginIndx : =>  Output: first index 
; EndIndx : =>  Output: last index 
; NbTIMESTAMP : =>  Output: Number of TIMESTAMP object 
; PacketGroup : =>  Output: ROI group name
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOFREEINFOTIMESTAMPGRP,TIMESTAMPtype,TIMESTAMPname,BeginIndx,EndIndx
;+
; TIMESTAMPtype : =>  Input: TIMESTAMP format 
; TIMESTAMPname : =>  Input: TIMESTAMP name list 
; BeginIndx : =>  Input: first index 
; EndIndx : =>  Input: last index 

;-
FUNCTION PIOGETTIMESTAMPPKT_SIZE,MyGroup
;+
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCREATEMAPGRP,Groupname,Coordsys,Ordering,NSide
;+
; Groupname : =>  Input: name of the opened group
; Coordsys : =>  Input: Coordinate type ("CELESTIAL"|"GALACTIC"|"ECLIPTIC")
; Ordering : =>  Input: Ordering type ("RING"|"NESTED")
; NSide : =>  Input: Healpix nside

;-
FUNCTION PIOCREATEMAPOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENMAPGRP,object,mode
;+
; object : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITEMAPOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEMAPOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEMAPOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEMAPOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSEMAPGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOINFOMAPGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,MAPtype=MAPtype,MAPname=MAPname,NbMAP=NbMAP,ordering=ordering,coordsys=coordsys,Nside=Nside
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; MAPtype : =>  Output: MAP format 
; MAPname : =>  Output: MAP name list 
; NbMAP : =>  Output: Number of MAP object 
; ordering : =>  Output: RING or NESTED 
; coordsys : =>  Output: GALACTIC,  ECLIPTIC or CELESTIAL 
; Nside : =>  Output: resolution 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOORDERINGGRP,ordering,MyGroup
;+
; ordering : =>  Output: RING or NESTED 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCOORDSYSGRP,coordsys,MyGroup
;+
; coordsys : =>  Output: GALACTIC,  ECLIPTIC or CELESTIAL 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIONSIDEGRP,MyGroup
;+
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCREATEPBRGRP,Groupname,ROIname,RingSize,Init_Size
;+
; Groupname : =>  Input: name of the opened group
; ROIname : =>  Input: Ring Information group name
; RingSize : =>  Input: Ring size
; Init_Size : =>  Input: default size

;-
FUNCTION PIOCREATEPOIOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOCREATEPHIOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENPBRGRP,object,mode
;+
; object : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITEPOIOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPOIOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPOIOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPOIOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPHIOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPHIOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPHIOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEPHIOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSEPBRGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOOPENPOIOBJECT,MyObjectName,typein,Command,step,mode,MyGroup
;+
; MyObjectName : =>  Input: Object name 
; typein : =>  Input: output data type 
; Command : =>  Input: output data type 
; step : =>  Input: step of reading 
; mode : =>  Input: 'w' or 'r' mode 
; MyGroup : =>  Input: link with the opened group 

;-
FUNCTION PIOGETPOIOBJECT,data,MyObject
;+
; data : =>  Input: for FORTRAN and IDL wrapper only 
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOFLUSHPOIOBJECT,MyObject
;+
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOSETPOIOBJECT,Tab,MyObject
;+
; Tab : =>  Input: data table 
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOSEEKPOIOBJECT,SampleNumber,MyObject
;+
; SampleNumber : =>  Input: number of samples 
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOCLOSEPOIOBJECT,MyObject
;+
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOINFOPBRGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,POItype=POItype,POIname=POIname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbPOI=NbPOI,PHItype=PHItype,PHIname=PHIname,NbPHI=NbPHI,ROIGroup=ROIGroup,RingSize=RingSize
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; POItype : =>  Output: POI format 
; POIname : =>  Output: POI name list 
; BeginIndx : =>  Output: first index 
; EndIndx : =>  Output: last index 
; NbPOI : =>  Output: Number of POI object 
; PHItype : =>  Output: POI format 
; PHIname : =>  Output: POI name list 
; NbPHI : =>  Output: Number of POI object 
; ROIGroup : =>  Output: ROI group name
; RingSize : =>  Output: last index 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIORINGSIZEGRP,MyGroup
;+
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCREATEROIGRP,Groupname,Init_Size
;+
; Groupname : =>  Input: name of the group
; Init_Size : =>  Input: default size

;-
FUNCTION PIOCREATEROIOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENROIGRP,object,mode
;+
; object : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITEROIOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEROIOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEROIOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEROIOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSEROIGRP,object
;+
; object : =>  Input: The opened group pointer 

;-
FUNCTION PIOINFOROIGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,ROItype=ROItype,ROIname=ROIname,BeginIndx=BeginIndx,EndIndx=EndIndx,NbROI=NbROI
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; ROItype : =>  Output: ROI format 
; ROIname : =>  Output: ROI name list 
; BeginIndx : =>  Output: first index 
; EndIndx : =>  Output: last index 
; NbROI : =>  Output: Number of ROI object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCREATEVECTGRP,Groupname
;+
; Groupname : =>  Input: name of the opened group

;-
FUNCTION PIOCREATEVECTOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENVECTGRP,object,mode
;+
; object : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITEVECTOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEVECTOBJECTIDX,data,Index,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEVECTOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEVECTOBJECTIDXMASK,data,Index,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSEVECTGRP,MyGroup
;+
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOINFOVECTGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,VECTtype=VECTtype,VECTname=VECTname,NbVECT=NbVECT
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; VECTtype : =>  Output: VECT format 
; VECTname : =>  Output: VECT name list 
; NbVECT : =>  Output: Number of VECT object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOOPENVECTOBJECT,MyObjectName,typein,Command,step,mode,MyGroup
;+
; MyObjectName : =>  Input: Object name 
; typein : =>  Input: output data type 
; Command : =>  Input: output data type 
; step : =>  Input: step of reading 
; mode : =>  Input: openning mode ["r" or "w"] 
; MyGroup : =>  Input: link with the opened group 

;-
FUNCTION PIOCLOSEVECTOBJECT,MyObject
;+
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOGETVECTOBJECT,data,MyObject
;+
; data : =>  Input: for FORTRAN and IDL wrapper only 
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOFLUSHVECTOBJECT,MyObject
;+
; MyObject : =>  Input: object structure 

;-
FUNCTION PIOSETVECTOBJECT,TabIn,MyObject
;+
; TabIn : =>  Input: input data table 
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOSEEKVECTOBJECT,SampleNumber,MyObject
;+
; SampleNumber : =>  Input: Number of sample to skip 
; MyObject : =>  Input: Object structure 

;-
FUNCTION PIOCREATETAB2DGRP,Groupname,NumberofColumn
;+
; Groupname : =>  Input: name of the opened group
; NumberofColumn : =>  Input: Number of column

;-
FUNCTION PIOCREATETAB2DOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENTAB2DGRP,GroupName,mode
;+
; GroupName : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITETAB2DOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB2DOBJECTIDX,data,Index1,Index2,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index1 : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB2DOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB2DOBJECTIDXMASK,data,Index1,Index2,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index1 : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSETAB2DGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOINFOTAB2DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,TAB2Dtype=TAB2Dtype,TAB2Dname=TAB2Dname,NbTAB2D=NbTAB2D
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; TAB2Dtype : =>  Output: TAB format 
; TAB2Dname : =>  Output: TAB name list 
; NbTAB2D : =>  Output: Number of TAB object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOGETTAB2DLINES,FirstLine,LastLine,ObjectName,MyGroup
;+
; FirstLine : =>  Output: the first available line 
; LastLine : =>  Output: the last available line 
; ObjectName : =>  Input: the object name 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOGETTAB2DCOLUMNGRP,MyGroup
;+
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCREATETAB3DGRP,Groupname,NAXIS1,NAXIS2
;+
; Groupname : =>  Input: name of the opened group
; NAXIS1 : =>  Input: Number of column
; NAXIS2 : =>  Input: Number of column

;-
FUNCTION PIOCREATETAB3DOBJECT,Objectname,typein,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENTAB3DGRP,GroupName,mode
;+
; GroupName : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOCLOSETAB3DGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB3DOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB3DOBJECTIDX,data,Index,Index2,Index3,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; Index3 : =>  Input: Input index 3 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB3DOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITETAB3DOBJECTIDXMASK,data,Index,Index2,Index3,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; Index3 : =>  Input: Input index 3 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOGETTAB3DAXISGRP,NAXIS1,NAXIS2,MyGroup
;+
; NAXIS1 : =>  Output: get Axis 1 
; NAXIS2 : =>  Output: get Axis 2  
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOGETTAB3DDEEP,FirstAx,LastAx,ObjectName,MyGroup
;+
; FirstAx : =>  Output: get the first samples Axis 3 
; LastAx : =>  Output: get the last sample Axis 3 
; ObjectName : =>  Input: name of the object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOINFOTAB3DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,TAB3Dtype=TAB3Dtype,TAB3Dname=TAB3Dname,NbTAB3D=NbTAB3D
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; TAB3Dtype : =>  Output: TAB format 
; TAB3Dname : =>  Output: TAB name list 
; NbTAB3D : =>  Output: Number of TAB object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOGETCRVALOBJECT,CRValX,CRValY,Objectname,MyGroup
;+
; CRValX : =>  Output: Sky coord. reference point
; CRValY : =>  Output: Sky coord. reference point
; Objectname : =>  Input: Objectname 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOSETLONPOLEGRP,LonPole,MyGroup
;+
; LonPole : =>  Input: Radesys 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETLONPOLEGRP,LonPole,MyGroup
;+
; LonPole : =>  Output: Radesys 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOSETLATPOLEGRP,LatPole,MyGroup
;+
; LatPole : =>  Input: Radesys 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETLATPOLEGRP,LatPole,MyGroup
;+
; LatPole : =>  Output: Radesys 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOSETRADESYSGRP,RadSysN,MyGroup
;+
; RadSysN : =>  Input: Radesys 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETRADESYSGRP,RadSysN,MyGroup
;+
; RadSysN : =>  Output: Radesys 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOSETEQUINOXGRP,Equinox,MyGroup
;+
; Equinox : =>  Input: Equinox 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETEQUINOXGRP,Equinox,MyGroup
;+
; Equinox : =>  Output: Equinox 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOSETROTOBJECT,Rot[9],Objectname,MyGroup
;+
; Rot[9] : =>  Input: RotAng in degree 
; Objectname : =>  Input: Objectname 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETROTOBJECT,Rot[9],Objectname,MyGroup
;+
; Rot[9] : =>  Output: RotAng in degree 
; Objectname : =>  Input: Objectname 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOSETROTANGOBJECT,RotAng,Objectname,MyGroup
;+
; RotAng : =>  Input: RotAng in degree 
; Objectname : =>  Input: Objectname 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETROTANGOBJECT,RotAng,Objectname,MyGroup
;+
; RotAng : =>  Output: RotAng in degree 
; Objectname : =>  Input: Objectname 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETPROJGRP,Proj,MyGroup
;+
; Proj : =>  Output: Projection definition 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETCOORDSYSGRP,Coordsys,MyGroup
;+
; Coordsys : =>  Output: Coordinate system  
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETAXISGRP,NAXIS1,NAXIS2,NAXIS3,MyGroup
;+
; NAXIS1 : =>  Output: Axis size 
; NAXIS2 : =>  Output: Axis size 
; NAXIS3 : =>  Output: Axis size 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETPIXSIZEGRP,PixSize1,PixSize2,MyGroup
;+
; PixSize1 : =>  Output: Pixel size 
; PixSize2 : =>  Output: Pixel size 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETCRPIXGRP,CRPix1,CRPix2,MyGroup
;+
; CRPix1 : =>  Output: central pixel definition 
; CRPix2 : =>  Output: central pixel definition 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOGETPVGRP,PV,MyGroup
;+
; PV : =>  Output: Projection parameter value 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOCREATEIMG2DGRP,Groupname,Proj,Coordsys,NAXIS1,NAXIS2,PixSize1,PixSize2,CRPixX,CRPixY,PV
;+
; Groupname : =>  Input: name of the opened group
; Proj : =>  Input: Projection type 
; Coordsys : =>  Input: Coordsys type 
; NAXIS1 : =>  Input: Size of the first axis
; NAXIS2 : =>  Input: Size of the second axis
; PixSize1 : =>  Input: Pixel size 1 in degree
; PixSize2 : =>  Input: Pixel size 2 in degree
; CRPixX : =>  Input: Pixel coord. reference point
; CRPixY : =>  Input: Pixel coord. reference point
; PV : =>  Input: Optional parameter

;-
FUNCTION PIOCREATEIMG2DOBJECT,Objectname,typein,CRValX,CRValY,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; CRValX : =>  Input: Sky coordinate of the reference pixel 
; CRValY : =>  Input: Sky coordinate of the reference pixel 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOOPENIMG2DGRP,GroupName,mode
;+
; GroupName : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOWRITEIMG2DOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEIMG2DOBJECTIDX,data,Index,Index2,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEIMG2DOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEIMG2DOBJECTIDXMASK,data,Index,Index2,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSEIMG2DGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOINFOIMG2DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,IMG2Dtype=IMG2Dtype,IMG2Dname=IMG2Dname,NbIMG2D=NbIMG2D
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; IMG2Dtype : =>  Output: IMG format 
; IMG2Dname : =>  Output: IMG name list 
; NbIMG2D : =>  Output: Number of IMG object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOCREATEIMG3DGRP,GroupName,Proj,Coordsys,NAXIS1,NAXIS2,NAXIS3,PixSize1,PixSize2,CRPixX,CRPixY,PV
;+
; GroupName : =>  Input: name of the opened group
; Proj : =>  Input: Projection type 
; Coordsys : =>  Input: Coordsys type 
; NAXIS1 : =>  Input: Size of the first axis
; NAXIS2 : =>  Input: Size of the second axis
; NAXIS3 : =>  Input: Size of the third axis
; PixSize1 : =>  Input: Pixel size 1 in degree
; PixSize2 : =>  Input: Pixel size 2 in degree
; CRPixX : =>  Input: Pixel coord. reference point
; CRPixY : =>  Input: Pixel coord. reference point
; PV : =>  Input: Optional parameter

;-
FUNCTION PIOOPENIMG3DGRP,GroupName,mode
;+
; GroupName : =>  Input: name of the group 
; mode : =>  Input: opening mode ("r"|"w") 

;-
FUNCTION PIOCREATEIMG3DOBJECT,Objectname,typein,CRValX,CRValY,MyGroup
;+
; Objectname : =>  Input: name of the created object 
; typein : =>  Input: type of the created object 
; CRValX : =>  Input: Sky coordinate of the reference pixel 
; CRValY : =>  Input: Sky coordinate of the reference pixel 
; MyGroup : =>  Input: The opened group pointer   

;-
FUNCTION PIOWRITEIMG3DOBJECT,data,Objectname,type,command,MyGroup,INDEX=INDEX,MASK=MASK
;+
; data : =>  Input: Input table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEIMG3DOBJECTIDX,data,Index,Index2,Index3,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; Index3 : =>  Input: Input index 3 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEIMG3DOBJECTMASK,data,Mask,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Mask : =>  Input: Input Mask table 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOWRITEIMG3DOBJECTIDXMASK,data,Index,Index2,Index3,Mask,NbData,Objectname,type,command,MyGroup
;+
; data : =>  Input: Input table 
; Index : =>  Input: Input index 1 
; Index2 : =>  Input: Input index 2 
; Index3 : =>  Input: Input index 3 
; Mask : =>  Input: Input Mask table 
; NbData : =>  Input: Table and Index size 
; Objectname : =>  Input: name of the object 
; type : =>  Input: data type 
; command : =>  Input: write command 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOCLOSEIMG3DGRP,MyGroup
;+
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIOINFOIMG3DGRP,MyGroup,FLGname=FLGname,NbFLG=NbFLG,IMG3Dtype=IMG3Dtype,IMG3Dname=IMG3Dname,NbIMG3D=NbIMG3D
;+
; FLGname : =>  Output: FLG name list 
; NbFLG : =>  Output: Number of FLG object 
; IMG3Dtype : =>  Output: IMG format 
; IMG3Dname : =>  Output: IMG name list 
; NbIMG3D : =>  Output: Number of IMG object 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOOPENIMOFILE,File,openm
;+
; File : =>  Input:  File Name 
; openm : =>  Input:  access mode

;-
FUNCTION PIOCLOSEIMO,MyGroup
;+
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOINFOIMO,MajorVersion,MinorVersion,MyGroup
;+
; MajorVersion : =>  Output: Major version 
; MinorVersion : =>  Output: Minor version 
; MyGroup : =>  Input: the opened group 

;-
FUNCTION PIOINITTABLE2COL,NbPoint,X_Unit,Y_Unit,X_Name,Y_Name,Comment,X_Type,Y_Type,MyTable2col
;+
; NbPoint : =>  Input: Number of points 
; X_Unit : =>  Input: X Unit 
; Y_Unit : =>  Input: Y Unit 
; X_Name : =>  Input: X Name 
; Y_Name : =>  Input: Y Name 
; Comment : =>  Input: Comment 
; X_Type : =>  Input: Type 
; Y_Type : =>  Input: Type 
; MyTable2col : =>  Input: Link to the IMOTable2col structure

;-
FUNCTION PIODELETETABLE2COL,MyTable2col
;+
; MyTable2col : =>  Input: Link to the IMOTable2col structure

;-
FUNCTION PIOFILLTABLE2COL,X_Value,ERR_X_Value,Y_Value,ERR_Y_Value,Index,MyTable2col
;+
; X_Value : =>  Input: X value 
; ERR_X_Value : =>  Input: Error on X value 
; Y_Value : =>  Input: Y value 
; ERR_Y_Value : =>  Input: Error on Y value 
; Index : =>  Input: Index Value 
; MyTable2col : =>  Input: Link to the IMOTable2col structure

;-
FUNCTION PIOEXTRACTTABLE2COL,X_Value,ERR_X_Value,Y_Value,ERR_Y_Value,Type,Index,MyTable2col
;+
; X_Value : =>  Output: X value 
; ERR_X_Value : =>  Output: Error on X value 
; Y_Value : =>  Output: Y value 
; ERR_Y_Value : =>  Output: Error on Y value 
; Type : =>  Input: Output Type 
; Index : =>  Input: Index Value 
; MyTable2col : =>  Input: Link to the IMOTable2col structure

;-
FUNCTION PIOSETTABLE2COL,MyData,Name,X_Type,Y_Type,MyIMO
;+
; MyData : =>  Input: Table2col 
; Name : =>  Input: Table2col name 
; X_Type : =>  Output : Data output type 
; Y_Type : =>  Output : Data output type 
; MyIMO : =>  Input: the opened group  ( IMO, BOLO, Etc.) 

;-
FUNCTION PIOINITTABLE,NbPoint,Unit,Name,Comment,Typein,MyTable
;+
; NbPoint : =>  Input: Number of points 
; Unit : =>  Input: Unit 
; Name : =>  Input: Name 
; Comment : =>  Input: Comment 
; Typein : =>  Input: Type 
; MyTable : =>  Input: Link to the IMOTable structure

;-
FUNCTION PIODELETETABLE,MyTable
;+
; MyTable : =>  Input: Link to the IMOTable structure

;-
FUNCTION PIOFILLTABLE,Value,ERR_Value,Index,MyTable
;+
; Value : =>  Input: X value 
; ERR_Value : =>  Input: Error on X value 
; Index : =>  Input: Index Value 
; MyTable : =>  Input: Link to the IMOTable  structure

;-
FUNCTION PIOEXTRACTTABLE,Value,ERR_Value,Type,Index,MyTable
;+
; Value : =>  Output : X value 
; ERR_Value : =>  Output : Error on X value 
; Type : =>  Input: Output Type 
; Index : =>  Input: Index Value 
; MyTable : =>  Input: Link to the IMOTable structure

;-
FUNCTION PIOSETTABLE,MyData,Name,Typein,MyIMO
;+
; MyData : =>  Input: Table2col 
; Name : =>  Input: Table2col name 
; Typein : =>  Input: Data output type 
; MyIMO : =>  Input: the opened group ( IMO, BOLO, Etc.) 

;-
FUNCTION PIOGETSTRINGLIST,MyDataList,Name,Attribs,MyIMO
;+
; MyDataList : =>  Output : String 
; Name : =>  Input: String name 
; Attribs : =>  Input: String attributes
; MyIMO : =>  Input: the opened group  ( IMO, BOLO, Etc.) 

;-
FUNCTION PIOGETSTRING,MyData,Name,Attribs,MyIMO
;+
; MyData : =>  Output : String 
; Name : =>  Input: String name 
; Attribs : =>  Input: String attributes
; MyIMO : =>  Input: the opened group  ( IMO, BOLO, Etc.) 

;-
FUNCTION PIOSETSTRING,MyData,Name,Attribs,MyIMO
;+
; MyData : =>  Input: String 
; Name : =>  Input: String name 
; Attribs : =>  Input: String attributes
; MyIMO : =>  Input: the opened group ( IMO, BOLO, Etc.) 

;-
FUNCTION PIOOPENPARFILE,ParName
;+
; ParName : =>  Input: Paramaeter file name

;-
FUNCTION PIOCLOSEPARFILE,MyPar
;+
; MyPar : =>  Input: Link to file 

;-
FUNCTION PIOGETPARAM,Value,Type,ParamName,MyPar
;+
; Value : =>  Output: Returned value 
; Type : =>  Input: returned type 
; ParamName : =>  Input: Parameter name 
; MyPar : =>  Input: Link to file 

;-
FUNCTION PIOGETLISTSZ,ListName,MyPar
;+
; ListName : =>  Input: Parameter name 
; MyPar : =>  Input: Link to file 

;-
FUNCTION PIOGETLIST,Value,Type,ListName,Index,MyPar
;+
; Value : =>  Output: Returned value 
; Type : =>  Input: returned type 
; ListName : =>  Input: Parameter name 
; Index : =>  Input: List index 
; MyPar : =>  Input: Link to file 

;-
FUNCTION PIOGETBOOLEAN,ParamName,MyPar
;+
; ParamName : =>  Input: Parameter name 
; MyPar : =>  Input: Link to file 

;-
FUNCTION PIOREADKEYWORDGRP,Value,Comment,InKeyword,Type,MyGroup
;+
; Value : =>  Output: Keyword value 
; Comment : =>  Output: comment on keyword 
; InKeyword : =>  Input: Keyword name 
; Type : =>  Input: Keyword type 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOWRITEKEYWORDGRP,Value,Comment,InKeyword,type,MyGroup
;+
; Value : =>  Input: Keyword value 
; Comment : =>  Input: comment on keyword 
; InKeyword : =>  Input: Keyword name 
; type : =>  Input: Keyword type 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIODELETEKEYWORDGRP,InKeyword,MyGroup
;+
; InKeyword : =>  Input: Keyword name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOKEYWORDLISTGRP,Keyword,Type,MyGroup
;+
; Keyword : =>  Output: Keyword name table
; Type : =>  Output: Keyword type table 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOWRITEKEYWORDOBJECT,Value,Comment,InKeyword,Type,Object,MyGroup,INDEX=INDEX,MASK=MASK
;+
; Value : =>  Input: Keyword value 
; Comment : =>  Input: Keyword name 
; InKeyword : =>  Input: Keyword name 
; Type : =>  Input: Keyword type 
; Object : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOREADKEYWORDOBJECT,Value,Comment,InKeyword,Type,Object,MyGroup
;+
; Value : =>  Output: Keyword value 
; Comment : =>  Input: Keyword name 
; InKeyword : =>  Input: Keyword name 
; Type : =>  Input: Keyword type 
; Object : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOKEYWORDLISTOBJECT,Keyword,Type,Object,MyGroup
;+
; Keyword : =>  Output: Keyword name table
; Type : =>  Output: Keyword type table 
; Object : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOGETGRPNAME,GroupName,ObjectName
;+
; GroupName : =>  Output : Group name 
; ObjectName : =>  Input: Group name 

;-
FUNCTION PIOGETROINAME,RoiName,GroupName
;+
; RoiName : =>  Output: The roi name 
; GroupName : =>  Input: the group name 

;-
FUNCTION PIOCHECKGROUP,GroupName
;+
; GroupName : =>  Input: the group name 

;-
FUNCTION PIOCHECKOBJECT,ObjectName,MyGroup
;+
; ObjectName : =>  Input: Object name 
; MyGroup : =>  Input: opened group reference 

;-
FUNCTION PIOGETGROUPTYPE,Typein,GroupName
;+
; Typein : =>  Output: Object name 
; GroupName : =>  Input: group name 

;-
FUNCTION PIOCHECKTEMPORARY,ObjectName
;+
; ObjectName : =>  Input: Object name 

;-
FUNCTION PIOGETGRPLIST,GrpList,GrpType,Database
;+
; GrpList : =>  Output : Group name list 
; GrpType : =>  Output : Group type list 
; Database : =>  Input: Database name 

;-
FUNCTION PIOGETOBJECTLIST,ObjList,ObjType,MyGroup
;+
; ObjList : =>  Output : Object name list 
; ObjType : =>  Output : Object type list 
; MyGroup : =>  Input: Group pointer 

;-
FUNCTION PIOCHECKDATABASE,Database
;+
; Database : =>  Input: Database name 

;-
FUNCTION PIOGETDATABASENAME,Database,HTMLPath
;+
; Database : =>  Output : Database name 
; HTMLPath : =>  Input: HTML database path 

;-
FUNCTION PIOGETLINKBASE,DataBaseName,GroupName
;+
; DataBaseName : =>  Output : Database name 
; GroupName : =>  Input: Group name 

;-
FUNCTION PIODELETEOBJECT,Objectname,MyGroup
;+
; Objectname : =>  Input: name of the object 
; MyGroup : =>  Input: The opened group pointer 

;-
FUNCTION PIODELETEGROUP,GroupName
;+
; GroupName : =>  Input: name of the group 

;-
FUNCTION PIOINFOOBJECT,TOItype,Datatype,BeginIndx,EndIndx,Author,Date,MyObjectName,MyGroup
;+
; TOItype : =>  Output: TOI format 
; Datatype : =>  Output: Data format 
; BeginIndx : =>  Output: first index 
; EndIndx : =>  Output: last index 
; Author : =>  Output: Object creator name 
; Date : =>  Output: Creation date 
; MyObjectName : =>  Input: Object name 
; MyGroup : =>  Input: link with the opened group 

;-
FUNCTION PIODELETEREPOSITORY,RepositoryName
;+
; RepositoryName : =>  Input: name of the database 

;-
