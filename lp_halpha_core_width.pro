pro lp_halpha_core_width, sp, ll, corewidth=corewidth, midp_corewidth=midp_corewidth, show=show
;
; Measures core width of Halpha (absorption line)
; Follows Cauzzi et al. 2009A&A...503..577C : 
; half width of average of wing positions at +-0.9 A from line core
; (here: delta 0.9 A from what is defined as 0)
;
; sp : spectrum
; ll : wavelength positions relative to line center, in AA
; corewidth : width in AA (units ll)
; midp_corewidth : midpoint of corewidth (in AA), sort of measure of shift
;

 if keyword_set(show) eq 0 then show=0 
 if min(ll) gt -0.9 or max(ll) lt 0.9 then begin
     print, 'spectral range should be in AA and at least cover -0.9 to +0.9'
     print, ll
     stop
 endif 
 nl=n_elements(ll)
 lrange=max(ll) - min(ll) ; wavelength range
; Cauzzi et al. 
; average wing intensities at delta lambda = 900 mAA
 spmax=max(sp, maxp)  ; spectrum max & min
 spmin=min(sp, minp)
 corewidth=0. & midp_corewidth=0.
 if show then begin
     window, 1
     plot, ll, sp, /xstyle
 endif 
 if minp gt 1 and minp lt nl-2 then begin ;check sensible range for min
     m0=min(abs(ll +0.9), l0)      ;find closest wavelength positions to -0.9 A
     m1=min(abs(ll -0.9), l1)      ;find closest wavelength positions to +0.9 A
     if m0 lt 0.001 and m1 lt 0.001 then begin ;If suitably close to 0.9
         avgw=mean([sp[l0], sp[l1]])            ;mean intensity of wavelength positions
         if show then begin
             plots, ll[l0], sp[l0], psym=1
             plots, ll[l1], sp[l1], psym=1
         endif 
     endif else begin
     ; if +-0.9 is not covered in spectral range, interpolate
         sp0=interpol(sp[0:minp-1], ll[0:minp-1], -0.9)
         sp1=interpol(sp[minp+1:*], ll[minp+1:*], 0.9)
         avgw=mean([sp0,sp1])
         if show then begin
             plots, -0.9, sp0, psym=1
             plots, +0.9, sp1, psym=1
         endif 
     endelse 
     ll0=interpol(ll[0:minp-1], sp[0:minp-1], (avgw-spmin)/2.+spmin) ;finds wavelength position short-side
     ll1=interpol(ll[minp+1:*], sp[minp+1:*], (avgw-spmin)/2.+spmin) ;finds wavelength position short-side
     corewidth=ll1-ll0
     midp_corewidth=ll0+(ll1-ll0)/2.
     if finite(corewidth) eq 0 then begin
         corewidth=0.
         midp_corewidth=0.
     endif 
     if corewidth gt lrange then begin
         corewidth=lrange
         midp_corewidth=0.
     endif 
     if show then begin
         title='core width = '+strtrim(corewidth,2)+' midpoint = '+strtrim(midp_corewidth,2)
         xyouts, (!x.window[1]-!x.window[0])/2. + !x.window[0], !y.window[1]+.01, title, alignment=0.5, /normal
         plots, ll0, (avgw-spmin)/2.+spmin, psym=1
          plots, ll1, (avgw-spmin)/2.+spmin, psym=1
          plots, [ll0, ll1], [1,1]+(avgw-spmin)/2.+spmin, li=1
          plots, ll[minp], spmin, psym=1
          plots, ll[maxp], spmax, psym=1
          plots, [1,1]*midp_corewidth, [!y.crange[0], (avgw-spmin)/2.+spmin], li=1
      endif 
  endif 
 return
end
