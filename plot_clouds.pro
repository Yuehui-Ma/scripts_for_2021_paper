pro plot_clouds 
catalog = '~/G105_overview/step7/clouds_12CO/ppb50/clusters_cen.cat'
clusters = '~/G105_overview/step7/clouds_12CO/ppb50/clusters.fits'
m0file = '~/G105_overview/step5/loc_U_m0.fits'
m0picfile = '~/G105_overview/step7/clouds_12CO/layer1_m0.eps'
lvfile = '~/G105_overview/step7/smoothcube/G105U_reb_clip_lvm0.fits'
lvpicfile = '~/G105_overview/step7/clouds_12CO/layer1_lv.eps';'~/My_work/G105_overview/step10/clouds_pv.eps'
vrange = [-110, -60]

cgloadct, 49
stretch, 0, 255, 0.8

fits_read, m0file, m0, m0hdr
m0 = sqrt(m0)
max = max(m0)
min = min(m0[m0 gt 0])*0.5
print, min, max
s = size(m0, /dim)
pic_aspect_ratio = double(s[0])/s[1] 
pos = [0.15, 0.15, 0.85, 0.85]
fits_read, clusters, data, hdr
nx = SXPAR(hdr,'NAXIS1')
ny = SXPAR(hdr,'NAXIS2')
xyad, hdr, [0, 0, nx-1, nx-1], [0, ny-1, ny-1, 0], l, b 
x_l = max(l)
x_u = min(l)
y_l = min(b)
y_u = max(b)

cgps_open, m0picfile, xsize = 9*pic_aspect_ratio, ysize = 9, /encapsulated, /portrait  
!p.thick = 2
!p.charthick = 1
!p.CHARSIZE = 2
img = bytscl(m0, min=min, max=max*0.8)
cgimage, img, position=pos, /save, /KEEP_ASPECT_RATIO
cgplot, 0, 0, /nodata, xrange =[x_l, x_u], yrange=[y_l, y_u], position=pos, /noerase,$
       xtitle = 'Galactic Longitude (!Uo!N)', ytitle = 'Galactic Latitude (!Uo!N)'
cgcolorbar, position = [pos[2]+0.015, pos[1], pos[2]+0.035, pos[3]], range = [min, max*0.8], /vertical, /right, $
            title = '(K km s!U-1!N)!U1/2!N'

readcol, catalog, id, v, l, b, format = 'd, d, d, d'
id = id[where(v/1000d gt 20)]
device, decompose=1
c_color = dblarr(n_elements(id))
if id[0] ne 0 then begin  
for i=0, n_elements(id)-1 do begin 
msk = dblarr(nx, ny)
loc = array_indices(data, where(data eq id[i]))
loc_x = loc[0, *]
loc_y = loc[1, *]
    for j=0, n_elements(loc_x)-1 do begin 
    msk[loc_x[j], loc_y[j]] = id[i]
    endfor
a = fix(randomu(seed)*255)
b = fix(randomu(seed)*255)
c = fix(randomu(seed)*255)
c_color[i] =  cgcolor([a, b ,c]) 
level = id[i]-1
cgcontour, msk, position = pos, /onimage, levels =level, label=0,  c_color = c_color[i]
endfor
endif
device, /close
cgps_close

;--------------------------------pv--------------------------------------------------
if id[0] ne 0 then begin  
device, decompose = 0

cgloadct, 49;, clip=[0,240]
stretch, 0, 255, 0.8

fits_read, lvfile , lv, m0hdr
fits_read, clusters, data, hdr
nx = SXPAR(hdr,'NAXIS1')
ny = SXPAR(hdr,'NAXIS3')
l = (dindgen(sxpar(hdr, 'NAXIS1')) + 1 - sxpar(hdr, 'CRPIX1')) * sxpar(hdr, 'CDELT1') + sxpar(hdr, 'CRVAL1')
v = (dindgen(sxpar(hdr, 'NAXIS3')) + 1 - sxpar(hdr, 'CRPIX3')) * sxpar(hdr, 'CDELT3') + sxpar(hdr, 'CRVAL3')
yr = (vrange - sxpar(hdr, 'CRVAL3')/1000.)/(sxpar(hdr, 'CDELT3')/1000.) + sxpar(hdr, 'CRPIX3') - 1
v = v[yr[0]:yr[1]]
data = data[*, *, yr[0]:yr[1]]
if (abs(sxpar(hdr, 'CDELT3')) gt 100) then v = v/1000d
x_l = max(l)
x_u = min(l)
y_l = min(v)
y_u = max(v)
print, x_l, x_u, y_l, y_u
lv = lv[*, yr[0]: yr[1]]
lv = sqrt(lv)
max = max(lv)
min = min(lv[lv gt 0])*1.5
print, min, max
s = size(lv, /dim)
pic_aspect_ratio = s[0]/s[1] 
pos = [0.15, 0.15, 0.85, 0.85]

cgps_open, lvpicfile, xsize = pic_aspect_ratio, ysize = 8, /encapsulated, /portrait  
!p.thick = 2
!p.charthick = 1
!p.CHARSIZE = 2
img = bytscl(smooth(lv, 2), min=min, max=max);/sqrt(8*alog(2))
cgimage, img, position=pos, /save;, /KEEP_ASPECT_RATIO
cgplot, 0, 0, /nodata, xrange =[x_l, x_u], yrange=[y_l, y_u], position=pos, /noerase,$
       xtitle = 'Galactic Longitude (!Uo!N)', ytitle = 'Velocity (km s!U-1!N)'
cgcolorbar, position = [pos[2]+0.015, pos[1], pos[2]+0.035, pos[3]], range = [min, max*0.8], /vertical, /right, $
            title = '(K arcdeg)!U1/2!N'
device, decompose=1
for i=0, n_elements(id)-1 do begin 
msk = dblarr(s[0], s[1])
loc = array_indices(data, where(data eq id[i]))
loc_x = loc[0, *]
loc_y = loc[2, *]
    for j=0, n_elements(loc_x)-1 do begin 
    msk[loc_x[j], loc_y[j]] = id[i]
    endfor
level = id[i]-1
cgcontour, msk, position = pos, /onimage, levels = level, label=0,  c_color = c_color[i]
endfor
device, /close
cgps_close
endif 
end 