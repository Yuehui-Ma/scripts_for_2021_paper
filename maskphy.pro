pro maskphy
name = ['loc', 'per', 'out']
path = '../step5/'
for i = 0, n_elements(name)-1 do begin 
    ;-----------------excitation temperature----------------------------
    fits_read, path+name[i]+'_U_peak.fits', peak12, hdr12
    T0 = 5.53213817d
    JTbg = (exp(T0/2.7)-1)^(-1)
    Tex = T0*(alog(1 + (peak12/T0 + JTbg)^(-1)))^(-1)
    Tex[where(Tex eq 2.7)] = 0
    fits_write, name[i]+'_Tex.fits', Tex, hdr12
    ;-----------------optical depth of 13CO------------------------------
    fits_read, path+name[i]+'_L_peak.fits', peak13, hdr13
    J13 = (exp(5.29d/Tex) - 1)^(-1)
    tau13 = -alog(1 - peak13/(5.29*(J13 - 0.164)))
    tau13[where(peak13 eq 0 or tau13 lt 0)] = 0
    tau13[where(Tex eq 0)] = 0
    fits_write, name[i]+'L_Tau.fits', tau13, hdr13
    ;-----------------optical depth of C18O------------------------------
    fits_read, path+name[i]+'_L2_peak.fits', peak18, hdr18
    J18 = (exp(5.27d/Tex) - 1)^(-1)
    tau18 = -alog(1 - peak18/(5.27*(J18 - 0.166)))
    tau18[where(peak18 eq 0 or tau18 lt 0)] = 0
    tau18[where(Tex eq 0)] = 0
    tau18[where(tau13 eq 0)] = 0
    fits_write, name[i]+'L2_Tau.fits', tau18, hdr18
    ;------------------column density------------------------------------
    ;---XCO---
    fits_read, path+name[i]+'_U_m0.fits', m012, hdr
    N12 = 2e20*m012/(1.1e4)
    fits_write, name[i]+'_N12.fits', N12, hdr
    ;---13CO---
    fits_read, path+name[i]+'_L_m0.fits', m013, hdr
    N13 = 2.42d14*m013/(1-exp(-5.29d/tex))*(1 + 0.88/tex) * tau13/(1-exp(-tau13))
    N13[where(finite(N13, /nan))] = 0
    N13[where(Tex eq 0)] = 0
    fits_write, name[i]+'_N13.fits', N13, hdr
    ;---C18O---
    fits_read, path+name[i]+'_L2_m0.fits', m018, hdr
    N18 = 2.54d14*m018/(1-exp(-5.27d/tex))*(1 + 0.88/tex) * tau18/(1-exp(-tau18))
    N18[where(finite(N18, /nan))] = 0
    ;N18[where(tau18 eq 0)] = 0
    N18[where(N13 eq 0)] = 0
    N18[where(N12 eq 0)] = 0
    fits_write, name[i]+'_N18.fits', N18, hdr

    openw, lun, name[i]+'_info.dat', /get_lun
    printf, lun, 'parameter  ', 'min  ','max  ', 'median  ', 'mean  ' 
    printf, lun, 'Tex ', min(Tex[where(Tex gt 0)]), max(Tex[where(Tex gt 0)]), $
            median(Tex[where(Tex gt 0)]), mean(Tex[where(Tex gt 0)]), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    printf, lun, 'Tex2 ', min(Tex[where(N13 gt 0)]), max(Tex[where(N13 gt 0)]), $
            median(Tex[where(N13 gt 0)]), mean(Tex[where(N13 gt 0)]), format = '(a, f11.2, f11.2, f11.2, f11.2)'            
    printf, lun, 'Tau13 ', min(Tau13[where(Tau13 gt 0)]), max(Tau13[where(Tau13 gt 0)]), $
            median(Tau13[where(Tau13 gt 0)]), mean(Tau13[where(Tau13 gt 0)]), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    if (where(Tau18 gt 0))[0] ne -1 then begin 
    printf, lun, 'Tex3 ', min(Tex[where(N18 gt 0)]), max(Tex[where(N18 gt 0)]), $
            median(Tex[where(N18 gt 0)]), mean(Tex[where(N18 gt 0)]), format = '(a, f11.2, f11.2, f11.2, f11.2)'    
    printf, lun, 'Tau13_2 ', min(Tau13[where(N18 gt 0)]), max(Tau13[where(N18 gt 0)]), $
            median(Tau13[where(N18 gt 0)]), mean(Tau13[where(N18 gt 0)]), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    printf, lun, 'Tau18 ', min(Tau18[where(Tau18 gt 0)]), max(Tau18[where(Tau18 gt 0)]), $
            median(Tau18[where(Tau18 gt 0)]), mean(Tau18[where(Tau18 gt 0)]), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    endif else begin 
    printf, lun, 'There is no effective C18O emission in the' + name[i] + ' spiral arm.'
    endelse 
    printf, lun, 'N12 ', min(alog10(N12[where(N12 gt 0)])), max(alog10(N12[where(N12 gt 0)])), $
            median(alog10(N12[where(N12 gt 0)])), mean(alog10(N12[where(N12 gt 0)])), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    printf, lun, 'N12_2 ', min(alog10(N12[where(N13 gt 0)])), max(alog10(N12[where(N13 gt 0)])), $
            median(alog10(N12[where(N13 gt 0)])), mean(alog10(N12[where(N13 gt 0)])), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    printf, lun, 'N13 ', min(alog10(N13[where(N13 gt 0)])), max(alog10(N13[where(N13 gt 0)])), $
            median(alog10(N13[where(N13 gt 0)])), mean(alog10(N13[where(N13 gt 0)])), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    if (where(N18 gt 0))[0] ne -1 then begin 
    printf, lun, 'N12_3 ', min(alog10(N12[where(N18 gt 0)])), max(alog10(N12[where(N18 gt 0)])), $
            median(alog10(N12[where(N18 gt 0)])), mean(alog10(N12[where(N18 gt 0)])), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    printf, lun, 'N13_3 ', min(alog10(N13[where(N18 gt 0)])), max(alog10(N13[where(N18 gt 0)])), $
            median(alog10(N13[where(N18 gt 0)])), mean(alog10(N13[where(N18 gt 0)])), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    printf, lun, 'N18 ', min(alog10(N18[where(N18 gt 0)])), max(alog10(N18[where(N18 gt 0)])), $
            median(alog10(N18[where(N18 gt 0)])), mean(alog10(N18[where(N18 gt 0)])), format = '(a, f11.2, f11.2, f11.2, f11.2)'
    endif else begin 
    printf, lun, 'There is no effective C18O emission in the' + name[i] + ' spiral arm.'
    endelse 
    free_lun, lun 
endfor 
    
cgps_open, 'MaskTex_stat.eps', /portrait, xsize = 8, ysize = 6, /encapsulated
!p.CHARSIZE = 2
!P.charthick = 2
!P.Thick = 3
!X.Thick = 3
!Y.Thick = 3
!Z.Thick = 3
pos = [0.2, 0.2, 0.9, 0.9]
color = ['red5', 'grn5', 'blu7']
label = ['Local', 'Perseus', 'Outer+OSC']    
cgplot, 0., 0, /nodata, xrange = [0, 40], yrange = [0, 0.3], xtitle = 'Tex (K)', ytitle = 'Relative Frequency', $
   xstyle = 1, ystyle = 1, position = pos 
;------------------------Tex-------------------------------------------
for i = 0, n_elements(name)-1 do begin 
    fits_read, name[i]+'_Tex.fits', Tex, hdr12 
    ;fits_read, name[i]+'_N13.fits', N13;, hdr12 
    ;help, Tex[where(Tex gt 0)]
    cghistoplot, Tex[where(Tex gt 0)], binsize = 0.5, position =  pos[*, 0], /oplot, $
       datacolorname = color[i], polycolor = color[i], thick = 5, $
       orientation = 45, /outline, /line_fill, /frequency
endfor 
items = label
sym = 27
al_legend, items, psym=sym, colors = color, /top, /right, charsize = 2, $
          charthick = 4, thick = 4, box = 0, symsize = 2
cgps_close

;------------------------Tau-------------------------------------------
cgps_open, 'MaskTau_stat.eps', /portrait, xsize = 8, ysize = 6, /encapsulated
!p.CHARSIZE = 2
!P.charthick = 2
!P.Thick = 3
!X.Thick = 3
!Y.Thick = 3
!Z.Thick = 3
pos = [0.2, 0.2, 0.9, 0.9]
color = ['red5', 'grn5']
label = ['Local', 'Perseus']        
cgplot, 0., 0, /nodata, xrange = [0, 2], yrange = [0, 0.15], xtitle = textoidl('\tau_{13CO}'), ytitle = 'Relative Frequency', $
   xstyle = 1, ystyle = 1, position =  pos
;------------------------Tau13-------------------------------------------
for i = 0, n_elements(name)-2 do begin 
    fits_read, name[i]+'L_Tau.fits', Tau13, hdr
    cghistoplot, Tau13[where(Tau13 gt 0)], binsize = 0.03, position =  pos[*, 0], /oplot, $
       datacolorname = color[i], polycolor = color[i], thick = 5, $
       orientation = 45, /outline, /line_fill, /frequency
endfor 
items = label
sym = 27
al_legend, items, psym=sym, colors = color, /top, /right, charsize = 2, $
          charthick = 4, thick = 4, box = 0, symsize = 2
cgps_close

cgps_open, 'MaskTau18_stat.eps', /portrait, xsize = 8, ysize = 6, /encapsulated
!p.CHARSIZE = 2
!P.charthick = 2
!P.Thick = 3
!X.Thick = 3
!Y.Thick = 3
!Z.Thick = 3
pos = [0.2, 0.2, 0.9, 0.9]
color = ['red5', 'grn5']
label = ['Local', 'Perseus']       
cgplot, 0., 0, /nodata, xrange = [0, 1], yrange = [0, 0.4], xtitle = textoidl('\tau_{C18O}'), ytitle = 'Relative Frequency', $
   xstyle = 1, ystyle = 1, position = pos
;------------------------Tau18-------------------------------------------
for i = 0, n_elements(name)-1 do begin 
    fits_read, name[i]+'L2_Tau.fits', Tau18, hdr
    if (where(Tau18 gt 0))[0] ne -1 then begin 
    cghistoplot, Tau18[where(Tau18 gt 0)], binsize = 0.03, position =  pos[*, 0], /oplot, $
       datacolorname = color[i], polycolor = color[i], thick = 5, $
       orientation = 45, /outline, /line_fill, /frequency
    endif
endfor 
items = label
sym = 27
al_legend, items, psym=sym, colors = color, /top, /right, charsize = 2, $
          charthick = 4, thick = 4, box = 0, symsize = 2
cgps_close


;------------------column density 12CO------------------------------
cgps_open, 'MaskN12_stat.eps', /portrait, xsize = 8, ysize = 6, /encapsulated
!p.CHARSIZE = 2
!P.charthick = 2
!P.Thick = 3
!X.Thick = 3
!Y.Thick = 3
!Z.Thick = 3
pos = [0.2, 0.2, 0.9, 0.9]
color = ['red5', 'grn5', 'blu7']
label = ['Local', 'Perseus', 'Outer+OSC']     
cgplot, 0., 0, /nodata, xrange = [15.5, 19], yrange = [0, 0.08], xtitle = textoidl('lg N_{12CO}'), ytitle = 'Relative Frequency', $
   xstyle = 1, ystyle = 1, position =  pos
;------------------------N12-------------------------------------------
for i = 0, n_elements(name)-1 do begin 
    fits_read, name[i]+'_N12.fits', N12, hdr 
    cghistoplot, alog10(N12[where(N12 gt 0)]), binsize = 0.05, position =  pos[*, 0], /oplot, $
       datacolorname = color[i], polycolor = color[i], thick = 5, /frequency, /line_fill, orientation = 40, /outline
endfor 
items = label
sym = 27
al_legend, items, psym=sym, colors = color, /top, /right, charsize = 2, $
          charthick = 4, thick = 4, box = 0, symsize = 2
cgps_close

;------------------column density 13CO------------------------------
cgps_open, 'MaskN13_stat.eps', /portrait, xsize = 8, ysize = 6, /encapsulated
!p.CHARSIZE = 2
!P.charthick = 2
!P.Thick = 3
!X.Thick = 3
!Y.Thick = 3
!Z.Thick = 3
pos = [0.2, 0.2, 0.9, 0.9]
color = ['red5', 'grn5']
label = ['Local', 'Perseus']     
cgplot, 0., 0, /nodata, xrange = [13.5, 18], yrange = [0, 0.12], xtitle = textoidl('lg N_{13CO}'), ytitle = 'Relative Frequency', $
   xstyle = 1, ystyle = 1, position =  pos
;------------------------N12-------------------------------------------
for i = 0, n_elements(name)-2 do begin 
    fits_read, name[i]+'_N13.fits', N13, hdr
    cghistoplot, alog10(N13[where(N13 gt 0)]), binsize = 0.1, position =  pos[*, 0], /oplot, $
       datacolorname = color[i], polycolor = color[i], thick = 5, /frequency, /line_fill, orientation = 40, /outline
endfor 
items = label
sym = 27
al_legend, items, psym=sym, colors = color, /top, /right, charsize = 2, $
          charthick = 4, thick = 4, box = 0, symsize = 2
cgps_close

;--------------------column density C18O--------------------------------------
cgps_open, 'MaskN18_stat.eps', /portrait, xsize = 8, ysize = 6, /encapsulated
!p.CHARSIZE = 2
!P.charthick = 2
!P.Thick = 3
!X.Thick = 3
!Y.Thick = 3
!Z.Thick = 3
pos = [0.2, 0.2, 0.9, 0.9]
color = ['red5', 'grn5']
label = ['Local', 'Perseus']    
cgplot, 0., 0, /nodata, xrange = [14, 16.4], yrange = [0, 0.15], xtitle = textoidl('lg N_{C18O}'), ytitle = 'Relative Frequency', $
   xstyle = 1, ystyle = 1, position =  pos
;------------------------N12-------------------------------------------
for i = 0, n_elements(name)-1 do begin 
    fits_read, name[i]+'_N18.fits', N18, hdr 
    if (where(N18 gt 0))[0] ne -1 then begin 
    cghistoplot, alog10(N18[where(N18 gt 0)]), binsize = 0.08, position =  pos[*, 0], /oplot, $
       datacolorname = color[i], polycolor = color[i], thick = 5, /frequency, /line_fill, orientation = 40, /outline
    endif 
endfor 
items = label
sym = 27
al_legend, items, psym=sym, colors = color, /top, /right, charsize = 2, $
          charthick = 4, thick = 4, box = 0, symsize = 2
cgps_close
end 
