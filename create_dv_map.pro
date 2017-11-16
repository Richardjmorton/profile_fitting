;PURPOSE - Fitting of line profiles using a variety of methods
;          Main purpose is the measurement of Doppler Velocity but
;          does provide profile minimum intensities
;          Wrapper for fitting routine fm_vel_spt
; 
;INPUTS - time - string in crispex form, i.e., '09.38.20'
;         lambda - wavelength in string, i.e., '6563'   
;         inpath - string of directory in which file is stored 
;         _EXTRA - additional keywords to pass to fm_vel_spt 
;       - method - method of fitting the line profile
;                  1) Centroid method
;                  2) Polynomial Fits
;                  3) Analytic fit to parabola - 3 points about min (very fast)
;                  4) 3 term Gaussfit fit to 3 points about min
;                  5) 3 term mpfit fit to 5 points about min
;
;
;OPTIONAL INPUTS - outpath - string of directory to save outputs to.
;                  stokes - /stokes , lets the routine know it contains stokes vectors
;				 - split - how many cores to split between, default=1. Generates
;                          IDL child processes on available cores. Will be slower
;                          if only splitting over a few cores (~<4). Only avilable for method 4.
;
;
;OUTPUTS - saves the velocity cube and profile minimum cube to 
;          a directory. Default is inpath.
;CALLING - Method 3 is very quick and gives comparable looking velocities to 4/5
;          Method 4 can be quick if using the split keyword and you have a lot of
;          cores available on your local machine!
;
;          DO NOT USE 1 AND 5 AT THE MOMENT - NOT FULLY FUNCTIONING 
;
;TO DO - Create error catching
;
;
;HISTORY - Created by RJ Morton and K. Mooroogen
;
;


PRO create_dv_map,time=time,lambda=lambda,inpath=inpath,outpath=outpath,stokes=stokes,_EXTRA=extra,$
	              filename=filename

;inpath - should be of the form /path/to/file/cripsex.wav.hr.min.sec
;inpath='SST/'

IF n_elements(inpath) EQ 0 THEN BEGIN
    cd,current=inp
    inpath=inp+'/'
ENDIF
help,inp
IF n_elements(filename) GT 0 THEN file=filename $
ELSE BEGIN
  IF NOT keyword_set(stokes) THEN $
		file=inpath+'crispex.'+lambda+'.'+time+'.time_corrected.icube' $
		ELSE file=inpath+'crispex.stokes.'+lambda+'.'+time+'.time_corrected.icube'
ENDELSE
IF n_elements(outpath) EQ 0 THEN outpath=inpath

lp_header,file,header=hd,dims=dim,nx=nx,ny=ny,nt=nt

; Set to long type
nx*=1L & ny*=1L & nt*=1L

restore,inpath+'wav'+lambda+'.idl'
nw=n_elements(wav) ;number of wavelengths

;Use wavelengths to calculate # time frames 
;as nt from lp_header is time*wavelength*stokes
;check for stokes in header
IF keyword_set(stokes) THEN nt=nt/nw/4 ELSE nt=nt/nw

;Open data
openr,lun,file,/get_lun
dat=assoc(lun,intarr(nx,ny,nw,/nozer),512)

lamb=float(lambda) ;convert string to float
                           
;Loop for sending individual frames to profile fitting routine
FOR i=0,nt-1 DO BEGIN
        ;read in Stokes I from icube - much quicker than loading in cube to memory
        print,'Processing frame number '+strtrim(i,2)+' of '+ strtrim(nt,2)
        IF keyword_set(stokes) THEN in=reform(dat[*,*,*,4.*i])*1. ELSE in=reform(dat[*,*,*,i])*1.
         
        ;Send the cube(x,y,lambda) for fitting  
		fm_vel_spt,in,wav,lambda=lamb,vel=a,intc=b,_EXTRA=extra

		;store
		IF n_elements(vel_map) EQ 0 THEN vel_map=a ELSE vel_map=[[[temporary(vel_map)]],[[a]]]
        IF n_elements(prof_min) EQ 0 THEN prof_min=b ELSE prof_min=[[[temporary(prof_min)]],[[b]]]
ENDFOR

;stop

;save velocity and profile minimum
save,vel_map,filename=outpath+'vel_map.sav'
IF n_elements(prof_min) GT 0 THEN save,prof_min,filename=outpath+'prof_min.sav'
END
