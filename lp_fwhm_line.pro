pro lp_fwhm_line, sp, ll, fwhm=fwhm, midp_fwhm=midp_fwhm, show=show
;
; Measure Full Width Half Maximum of an absorption line
; maximum : range between minimum and maximum in spectrum
; Note: not adapted to deal with emission line
;
; sp : spectrum
; ll : wavelength positions
; fwhm : FWHM in units of ll
; midp_fwhm : midpoint in wavelength units of FWHM, sort of measure of shift
; show : 1 : show results in plot
;

 if keyword_set(show) eq 0 then show=0 
 nl=n_elements(ll)
 lrange=max(ll) - min(ll) ; wavelength range
 spmax=max(sp, maxp)      ; spectrum max & min
 spmin=min(sp, minp)
 sprange=spmax-spmin      ;spectrum range
 fwhm=0. & midp_fwhm=0.
 if show then begin
     window, 1
     plot, ll, sp, /xstyle 
 endif 
 if minp gt 1 and minp lt nl-2 then begin
     ll0=interpol(ll[0:minp-1], sp[0:minp-1], sprange/2.+spmin) ;finds wavelength position at half maximum short-side
     ll1=interpol(ll[minp+1:*], sp[minp+1:*], sprange/2.+spmin) ;finds wavelength position at half maximum long-side
     fwhm=ll1-ll0
     midp_fwhm=ll0+(ll1-ll0)/2.
     if finite(fwhm) eq 0 then begin
         fwhm=0.
         midp_fwhm=0.
     endif 
     if fwhm gt lrange then begin  ; FWHM cannot be larger than wavelength coverage
         fwhm=lrange
         midp_fwhm=0.
     endif 
     if show then begin
          title='FWHM = '+strtrim(fwhm,2)+' midpoint = '+strtrim(midp_fwhm,2)
          xyouts, (!x.window[1]-!x.window[0])/2. + !x.window[0], !y.window[1]+.01, title, alignment=0.5, /normal
          plots, ll0, sprange/2.+spmin, psym=1
          plots, ll1, sprange/2.+spmin, psym=1
          plots, [ll0, ll1], [1,1]+sprange/2.+spmin
          plots, ll[minp], spmin, psym=1
          plots, ll[maxp], spmax, psym=1
          plots, [1,1]*midp_fwhm, [!y.crange[0], sprange/2.+spmin], li=1
      endif
  endif 
;stop
 return
end
