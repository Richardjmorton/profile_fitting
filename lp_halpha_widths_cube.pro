;INPUTS:
;   imcubename - path to image cube
;   spcubename - path to sp cube
;   nl         - number of wavelengths (required if no spectfile) 
;   nt 	     - number of time-steps  (required if no sp cube) 
;   spectfilename - path to file containing offset wave-lengths (idl save file)
;   ll          - wavelength offsets (required if no spectfile) 
;   pref       - reference wavelength (default halpha)
;   outcubename - path to save results to (will use imcubename if not specified)
;   fwhm_only - calculate width at FWHM
;   corewidth_only - calculate width for line core only
;show=show


pro lp_halpha_widths_cube, imcubename=imcubename, spcubename=spcubename, nl=nl, nt=nt, $
  spectfilename=spectfilename, ll=ll, pref=pref, outcubename=outcubename, $
  fwhm_only=fwhm_only, corewidth_only=corewidth_only, show=show

 if keyword_set(fwhm_only) eq 0 then fwhm_only=0
 if keyword_set(corewidth_only) eq 0 then corewidth_only=0
 if keyword_set(pref) eq 0 then pref=6563  ; default is Halpha
 if keyword_set(show) eq 0 then show=0

 if keyword_set(spcubename) eq 0 then begin
     if keyword_set(nt) eq 0 then begin
         print, 'provide number of time steps, nt='
         stop
     endif
 endif else lp_header, spcubename, nx=nl, ny=nt
 
 if keyword_set(spectfilename) eq 0 then begin
     if keyword_set(ll) eq 0 then begin
         print, 'provide offset wavelengths ll='
         stop
     endif 
 endif else begin
     restore, spectfilename, /verbose
     if keyword_set(nl) eq 0 then nl=n_elements(spect_pos)
     m=min(abs(spect_pos - float(pref)), lc)
     print, 'line center at '+strtrim(lc,2)+' spect_pos = '+strtrim(spect_pos[lc],2)
     ll=spect_pos - float(pref)
 endelse 
 
 ;Extracts data parameters from header
 lp_header, imcubename, nx=nx, ny=ny
 if keyword_set(outcubename) eq 0 then begin
     p=strpos(imcubename, '/', /reverse_search)
     dir=strmid(imcubename, 0, p)+'/'
     t=strsplit(imcubename, '/', /extract)
     fn=t[-1]
     p=strpos(fn, '.')
     if fwhm_only eq 0 and corewidth_only eq 0 then outcubename=dir+'widths'+strmid(fn, p)
     if fwhm_only then outcubename=dir+'fwhm'+strmid(fn, p)
     if corewidth_only then outcubename=dir+'corewidth'+strmid(fn, p)
 endif 
; results shoud be integer -> .icube
 t=strsplit(outcubename, '.', /extract)
 ext=t[-1]
 if ext ne 'icube' then begin
     p=strpos(outcubename, '.', /reverse_search)
     outcubename=strmid(outcubename, 0, p)+'.icube'
 endif 
 print, 'saving results in '+outcubename
stop
 if fwhm_only eq 0 and corewidth_only eq 0 then begin
     p=strpos(outcubename, '.', /reverse_search)
     infofilename=strmid(outcubename, 0, p)+'.info.txt'
     openw, luw, infofilename, /get_lun
     printf, luw, '0: FWHM [mAA]'
     printf, luw, '1: midpoint FWHM [mAA]'
     printf, luw, '0: core width [mAA]'
     printf, luw, '0: midpoint core width [mAA]'
     free_lun, luw
     print, 'wrote info in '+infofilename
 endif 

 openr, lur, imcubename, /get_lun
 r=assoc(lur, intarr(nx,ny,nl), 512)
 openw, luw, outcubename, /get_lun
 if fwhm_only eq 0 and corewidth_only eq 0 then h=make_lp_header(intarr(nx,ny), nt=nt*4) else h=make_lp_header(intarr(nx,ny), nt=nt)
 w=assoc(luw, bytarr(512)) & w[0]=byte(h)
 w=assoc(luw, intarr(nx,ny), 512)
 if show then window, 0, xs=nx,ys=ny, retain=2
 for t=0, nt-1 do begin
     cube=r[t]
     widthmap=intarr(nx,ny)
     shiftmap=intarr(nx,ny)
     cwidthmap=intarr(nx,ny)
     cshiftmap=intarr(nx,ny)
     for x=0, nx-1 do begin
         for y=0, ny-1 do begin
             sp=cube[x,y,*] 
             if corewidth_only eq 0 then begin
                 lp_fwhm_line, sp, ll, fwhm=fwhm, midp_fwhm=midp_fwhm, show=0
                 widthmap[x,y]=fix(round(fwhm*1000))
                 shiftmap[x,y]=fix(round(midp_fwhm*1000))
             endif 
             if fwhm_only eq 0 then begin
                 lp_halpha_core_width, sp, ll, corewidth=corewidth, midp_corewidth=midp_corewidth, show=0
                 cwidthmap[x,y]=fix(round(corewidth*1000))
                 cshiftmap[x,y]=fix(round(midp_corewidth*1000))
             endif 
         endfor 
         print,string(13b)+' % finished: ',float(x)*100./(nx-1),format='(a,f4.0,$)'
     endfor
     if fwhm_only eq 0 and corewidth_only eq 0 then begin
         w[t*4]=widthmap
         w[t*4+1]=shiftmap
         w[t*4+2]=cwidthmap
         w[t*4+3]=cshiftmap
     endif 
     if fwhm_only then w[t]=widthmap
     if corewidth_only then w[t]=cwidthmap
     if show then begin
         if fwhm_only eq 1 or corewidth_only eq 0 then tvscl, histo_opt(widthmap) else tvscl, histo_opt(cwidthmap)
     endif 
     print, strtrim(t,2)+'/'+strtrim(nt,2)
 endfor 
 free_lun, luw, lur
 print, 'wrote: '+outcubename
 if fwhm_only eq 0 and corewidth_only eq 0 then make_crispex_sp_cube, outcubename, 4, nt

return
end
