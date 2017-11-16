
;Centroid method
FUNCTION COM, ivg,wav_diff

  common constants, c, lambda ,nx,ny,nw 

  vel=total(ivg*wav_diff,3)/total(ivg,3)*c/lambda/1.e3

END

;Polynomial fit
FUNCTION polyf, ivg,wav_diff

  common constants, c, lambda,nx,ny,nw 

  res=min(ivg,dim=3,loc)
  vel=fltarr(nx,ny)
  x=transpose(rebin([-1,0,1],3,nx,ny),[1,2,0])*nx*ny+rebin(loc,nx,ny,3)
  y=ivg[x]
  wa=wav_diff[x]  

  FOR j=0,ny-1 DO FOR i=0,nx-1 DO BEGIN
       xin=reform(wa[i,j,*])
       yin=reform(y[i,j,*])
       
       res=poly_fit(xin,yin,2)
       vel[i,j]=-res[1]/2./res[2]   ;minimum of parabola
  ENDFOR

END

;Analytic minimum of parabola
FUNCTION parabmin,ivg,wav_diff

  common constants, c, lamb,nx,ny,nw 

  res=min(ivg,dim=3,loc)
  x=transpose(rebin([-1,0,1],3,nx,ny),[1,2,0])*nx*ny+rebin(loc,nx,ny,3)
  ;x2=transpose(rebin([-1,0,1],3,nx,ny),[1,2,0])+rebin(loc,nx,ny,3)/nx/ny
  y=ivg[x]
  wa=wav_diff[x]

  top=wa[*,*,2]^2*(y[*,*,1]-y[*,*,0])+wa[*,*,1]^2*(y[*,*,2]-y[*,*,0])+wa[*,*,0]^2*(y[*,*,1]-y[*,*,2])
  bot=wa[*,*,2]*(y[*,*,1]-y[*,*,0])+wa[*,*,1]*(y[*,*,2]-y[*,*,0])+wa[*,*,0]*(y[*,*,1]-y[*,*,2])

  vel=-0.5*top/bot
  in=where(loc/nx/ny gt nw-2)
  IF n_elements(in) gt 0 THEN vel[in]=wav_diff[loc[in]/nx/ny] 
  in=where(loc/nx/ny lt 1)
  IF n_elements(in) gt 0 THEN vel[in]=wav_diff[loc[in]/nx/ny]   

  return,vel
END

;Gaussian fit
FUNCTION gfit, ivg,wav_diff

  common constants, c, lamb,nx,ny,nw 

  res=min(ivg,dim=3,loc)
  x=transpose(rebin([-2,-1,0,1,2],5,nx,ny),[1,2,0])*nx*ny+rebin(loc,nx,ny,5)
  ;x2=transpose(rebin([-2,-1,0,1,2],5,nx,ny),[1,2,0])+rebin(loc,nx,ny,5)/nx/ny
  y=ivg[x]
  wa=wav_diff[x]

  vel=fltarr(nx,ny)
  err=fltarr(5)
  err[*]=1.

  FOR j=0,ny-1 DO FOR i=0,nx-1 DO BEGIN
      xin=reform(wa[i,j,*])
      yin=reform(y[i,j,*])
      yin-=max(yin)
      estimates=[min(yin),xin[2],0.1]
       

      res=gaussfit(xin,yin,nterms=3,coeff,estimates=estimates)
      vel[i,j]=coeff[1]  
      ;plot,xin,yin,psym=1
      ;oplot,xin,mygaussred(xin,coeff)
      ;stop
      ;counter,j*ny+i,nx*ny,/percent
  ENDFOR 

  ;Will definitely be bad values
  in=where(loc/nx/ny gt nw-2)
  IF n_elements(in) gt 0 THEN vel[in]=wav_diff[loc[in]/nx/ny] 
  in=where(loc/nx/ny lt 1)
  IF n_elements(in) gt 0 THEN vel[in]=wav_diff[loc[in]/nx/ny]   
  in=where(vel gt max(wav_diff))
  vel[in]=max(wav_diff)
  in=where(vel lt min(wav_diff))
  vel[in]=min(wav_diff)
  
  return,vel
END

;Gaussian fit
FUNCTION gfitm, ivg,wav_diff

  common constants, c, lamb,nx,ny,nw 

  res=min(ivg,dim=3,loc)
  x=transpose(rebin([-2,-1,0,1,2],5,nx,ny),[1,2,0])*nx*ny+rebin(loc,nx,ny,5)
  ;x2=transpose(rebin([-2,-1,0,1,2],5,nx,ny),[1,2,0])+rebin(loc,nx,ny,5)/nx/ny
  y=ivg[x]
  wa=wav_diff[x]

  vel=fltarr(nx,ny)
  err=fltarr(5)
  err[*]=1.

  FOR j=0,ny-1 DO FOR i=0,nx-1 DO BEGIN
      xin=reform(wa[i,j,*])
      yin=reform(y[i,j,*])
      yin-=max(yin)
      estimates=[min(yin),xin[2],0.1]
      weights=fltarr(n_elements(xin))
      weights[*]=1.       

      coeff=mpfitfun('mygaussred',xin,yin,weights,estimates,weights=weights,xtol=1e-3,/quiet)
      
      IF n_elements(coeff) gt 1 THEN vel[i,j]=coeff[1] 
      
      counter,j*ny+i,nx*ny,/percent
  ENDFOR 

  ;Will definitely be bad values
  in=where(loc/nx/ny gt nw-2)
  IF n_elements(in) gt 0 THEN vel[in]=wav_diff[loc[in]/nx/ny] 
  in=where(loc/nx/ny lt 1)
  IF n_elements(in) gt 0 THEN vel[in]=wav_diff[loc[in]/nx/ny]   
  in=where(vel gt max(wav_diff))
  vel[in]=max(wav_diff)
  in=where(vel lt min(wav_diff))
  vel[in]=min(wav_diff)
  
  return,vel
END

;COG method
FUNCTION cogfit, ivg,wav_diff

  common constants, c, lamb,nx,ny,nw 

  conto=smooth(ivg[*,*,0],9,/edge_truncate)
  contt=smooth(ivg[*,*,nw-1],9,/edge_truncate)
  icont=rebin((conto+contt)/2.,nx,ny,nw)
  wav_diff=wav_diff+lamb

  fint=wav_diff*(icont-ivg)
  sint=(icont-ivg)

  vel=fltarr(nx,ny)
  vel=total(fint,3)/total(sint,3)

  !p.multi=[0,2,2]
  ;FOR i=0,nx-1 DO FOR j=0,ny-1 DO BEGIN
      
  ;     po=int_tabulated(reform(wav_diff[i,j,0:nw-1]), reform( fint[i,j,*]))
  ;     pt=int_tabulated(reform( wav_diff[i,j,0:nw-1]), reform( sint[i,j,*]))
  ;     vel[i,j]=po/pt
     
       ;plot,ivg[i,j,*]
       ;plots,[0,15],[1,1]*icont[i,j,*],linestyle=2
       ;plot,fint[i,j,*]
       ;plot,sint[i,j,*]
       ;pause
  ;    counter,i*ny+j,nx*ny,/percent
  ;ENDFOR 

  ;Will definitely be bad values
  in=where(vel gt max(wav_diff))
  vel[in]=max(wav_diff)
  in=where(vel lt min(wav_diff))
  vel[in]=min(wav_diff)
  
  return,vel
END



;Find moments of spectral line to measure Doppler velocities
; cube - spectral cube
; wav - wavelengths
PRO fm_vel,cube,wav, lambda=lambda,method=method,vel=vel

;on_error,2 

common constants, c, lamb ,nx,ny,nw 

IF n_elements(lambda) EQ 0 THEN message,'Need central wavelength in A'

c=299792458*1d  ;speed of light
lamb=lambda

sz=size(cube)

IF sz(0) EQ 3 THEN BEGIN
   nx=sz(1) & ny=sz(2) & nw=sz(3)

ENDIF

;FOV average
avgprof=total(total(cube,1),1)/nx/ny
avgprof=transpose(rebin(temporary(avgprof),nw,nx,ny),[1,2,0])


wav_diff=wav ;(lambda-lambda_0)
wav_diff=transpose(rebin(temporary(wav_diff),nw,nx,ny),[1,2,0])

ivg=cube ;(I_avg-I)

CASE method OF 

    1: res=com(ivg,wav_diff)
    2: res=polyf(ivg,wav_diff)
    3: res=parabmin(ivg,wav_diff)
    4: res=gfit(ivg,wav_diff)
    5: res=cogfit(ivg,wav_diff)
    6: res=gfitm(ivg,wav_diff)

ENDCASE


vel=res/lamb*c


END
