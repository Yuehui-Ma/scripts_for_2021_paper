; Author : @Yuehui Ma
; Date : 2019-07
; For make panels of NPDFs from given H2 column density files, fit each N-PDF with LN and LN+PL functions using mpfit, and select a better model with a smaller reduced chi-squared value. 

pro pdf_panels
cgloadct,0
cat = '../step7/clouds_13CO/figures_g19/sort_clouds_area.dat'
fig = './figures/pdf_panel_loc.eps'
spawn, 'rm ' + fig
openw, lun, './figures/fit_parameters_loc.tex', /get_lun
readcol, cat, name, id, l, b, v, d, format = 'a, d, d, d, d, d'
sel = where(v gt -27)
binsize = 0.15   
xrange = [-4, 4]
yrange = [1e-4, 4]
ytickv = [1e-4, 1e-3, 1e-2, 1e-1, 1]
ytickname = ['10!U-4!N','10!U-3!N','10!U-2!N','10!U-1!N','10!U0!N']
cgps_open, fig,  xsize = 18, ysize = 16 , /encapsulated
!p.multi = [0, 5, 6]
!P.charsize = 2.3
!P.Thick = 3
!P.CharThick = 1
!X.Thick = 1
!Y.Thick = 1
!Z.Thick = 1
name = name[sel]
id = id[sel]
l = l[sel]
b = b[sel]
v = v[sel]
d = d[sel]
width = []
mach = []
for i = 0, n_elements(id) - 1 do begin 
    path = './'+name[i]+'/'
    file = path + name[i] + 'L_N.fits' 
    limfile = path + name[i]+'info.dat'
    phyinfo = path + name[i]+'_phyp.dat'
    readcol, phyinfo, nam, idd, Tkin, ms, radi, M, alpha, format = '(a, I, d, d, d, d, d)'
    fits_read, file, N
    good = where(N gt 0)
    dat = N[good]
    help, dat
    ;===================for table ============================
    N0 = mean(dat)
    s = alog(dat/mean(dat))
    mudata = mean(s)
    sdata = stdev(s)
    ;=========================================================
    ;print, min(dat), max(dat), mean(dat), median(dat)
    data = alog(dat/mean(dat))
    min = min(data)
    max = max(data)
    hist0 = [histogram(data, binsize = binsize, min = min, max = max, location = loc), 0]
    hist = float(hist0)/(float(total(hist0))*binsize)
    nhist = n_elements(hist)
    bins = lindgen(nhist)*binsize + min
    x = fltarr(2*nhist)
    x[2*lindgen(nhist)] = bins
    x[2*lindgen(nhist)+1] = bins 
    y = fltarr(2*nhist)
    y[2*lindgen(nhist)] = hist
    y[2*lindgen(nhist)+1] = hist
    y = shift(y, 1)
;------------------------------------------------------------------
    pos = [0.15, 0.2, 0.9, 0.9]
    plot, x, y, xTitle = 's = ln(N/<N>)', ytitle = 'p(s)', /ylog, $
      yminor = 9, xminor = 10, xstyle = 1+8, xrange = xrange, yrange = yrange,$
      xtick_get = tic, ystyle = 1, ytickv = ytickv, ytickname = ytickname, $
      xmargin = [8, 8], ytick_get = tic2, /nodata
    xmin = mean(dat)*exp(tic[0])
    xmax = mean(dat)*exp(tic[n_elements(tic)-1])
    axis, xaxis = 1, xstyle = 1, xrange = [xmin, xmax], /xlog
;----------gaussian fit-------------
    binCenters = bins[0:nhist-2] + (binsize / 2.0)
    mad = sqrt(float(hist0[0:nhist-2]))/(total(hist0[0:nhist-2])*binsize)
    cgplot, bincenters, hist[0:nhist-2], err_yhigh = mad, err_ylow = mad, err_color = 'indian red', $
        err_thick=1.5, err_width = 0.001, psym = 3, /overplot
    readcol, limfile, limit, format = 'd'
    mass = limit[1]
    cgplot, [alog(limit[2]/mean(dat)), alog(limit[2]/mean(dat))], [0.00001, 3], /overplot, thick = 2, linestyle = 3;, color = 'grey'
    numy = [10/(total(hist0)*binsize), 10/(total(hist0)*binsize)]
    cgplot, xrange, numy, thick = 2, linestyle = 1,  /overplot
    oplot, x, y, color = cgcolor('black')
    cgtext, -1.5, 0.3e-3, name[i], color = 'black', charsize = 1, charthick = 2
    range = bincenters[where(hist0 gt 10)]
    fitrange_g = [alog(min(dat)/mean(dat)), range[-1]] 
    n1 = where(bincenters ge fitrange_g[0] and bincenters le fitrange_g[1]);2.72
    dat1 = binCenters[n1]
    h0 = hist0[0:nhist-2]
    h = hist[0:nhist-2]
    dat2 = h[n1]
    measure_errors = sqrt(h0[n1])/(total(hist0)*binsize)
    parinfo = replicate({value:0.D, fixed:0, limited:[1, 1], limits:[0.D, 0.D]}, 3)
    parinfo[0].value = -0.5
    parinfo[1].value = 0.6
    parinfo[2].value = 1. 
    parinfo[0].limits = [-2., 2.]
    parinfo[1].limits = [0., 2.]
    parinfo[2].limits = [0.5, 2.]
    expr = 'GAUSS1(X, P)'
    start1 = [ -0.5, 0.6, 1.]
    result1 = MPFITEXPR(expr, dat1, dat2, measure_errors, perror = perror, quiet = 1, parinfo = parinfo, $
          bestnorm = bestnorm1, DOF = DOF1)
    PCERROR = PERROR * SQRT(BESTNORM1 / DOF1)
    red_chi1 = BESTNORM1 / DOF1
    
;--------------------------powerlaw fit----------------------------------
    binCenters = bins + (binsize / 2.0)
    left = where(hist eq max(hist))
    range = bincenters[where(hist0 gt 10)]
    fitrange_p = [bincenters[left], range[-1]]
    n1 = where(bincenters ge alog(4e21/mean(dat)) and bincenters le fitrange_p[1]);2.72
    dat1 = binCenters[n1]
    h0 = hist0
    h = alog(hist)
    dat2 = h[n1]
    measure_errors = (sqrt(h0[n1])/(total(hist0)*binsize))/hist[n1]
    parinfo = replicate({value:0.D, fixed:0, limited:[1, 1], limits:[0.D, 0.D]}, 2)
    parinfo[0].value = -2.0
    parinfo[1].value = -0.5
    parinfo[0].limits = [-3., -1.]
    parinfo[1].limits = [-1.5, 0.]
    expr = 'p[0]*X + p[1]'
    start2 = [-2.0, -0.5]
    result2 = MPFITEXPR(expr, dat1, dat2, measure_errors, perror = perror, quiet = 1, parinfo = parinfo, $
          bestnorm = bestnorm2, DOF = DOF2)
    PCERROR = PERROR * SQRT(BESTNORM2 / DOF2)
    red_chi2 = BESTNORM2 / DOF2
    if (name[i] eq 'G113.31-1.22') then continue
    if (finite(result1[0]) and ~finite(result2[0]) and (result1[0] ne start1[0]) and dof1 ne 0) then begin 
        ;=========================================for shadows====================================================
        cgpolygon, [fitrange_g[0], fitrange_g[1], fitrange_g[1], fitrange_g[0], fitrange_g[0]], $
          [yrange[0], yrange[0], yrange[1], yrange[1], yrange[0]], /data, color = 'dodger blue',linestyle = 1, thick=3;, /fill
        ;=======================================================================================================
        cgplot, bincenters, gauss1(bincenters, result1), color = 'dodger blue', /overplot, linestyle = 2, thick = 3
        center = String(result1[0], FORMAT='(F0.2)')
        sigma = String(result1[1], FORMAT='(F0.2)')
        text1 = textoidl('\mu = ')
        text2 = textoidl('\sigma = ')
        cgtext, 0.8, 0.6, /data, text1 + center, COLOR = 'blue', charsize = 1
        cgtext, 0.8, 0.2, /DATA, text2 + sigma, COLOR = 'blue', charsize = 1
        printf, lun, name[i], '&', d[i], '&', mass, '&', N0, '&', mudata, '&', sdata, '&', result1[0], $
               '&', result1[1], '&', '...', '&', 'LN', '&', '...', '\\', $
                format='(a,a,f7.2,a,e11.2,a,e11.2,a,f7.2,a,f7.2,a,f7.2,a,f7.2,a, a ,a, a ,a, a,a)' ;  log-normal fitted
        ;
        width = [width, result1[1]]
        mach = [mach, ms]
    endif 
    if (~finite(result1[0]) and finite(result2[0]) and (result2[0] ne start2[0]) and dof2 ne 0) then begin ;
        ;=========================================for shadows====================================================
        cgpolygon, [fitrange_p[0], fitrange_p[1], fitrange_p[1], fitrange_p[0], fitrange_p[0]], $
         [yrange[0], yrange[0], yrange[1], yrange[1], yrange[0]], /data, color = 'magenta',linestyle = 1, thick = 3;, /fill
        ;========================================================================================================= 
        cgplot, bincenters, exp(result2[0]*bincenters+result2[1]), color = 'magenta', /overplot, linestyle = 2, thick = 3
        sigma = String(result2[1], FORMAT='(F0.2)')
        text1 = textoidl('\mu = ')
        text2 = textoidl('\sigma = ')
        cgtext, 0.8, 0.2, /DATA, text2 + sigma, COLOR='magenta', charsize = 1
        printf, lun, name[i], '&', d[i], '&', mass, '&', N0, '&', mudata, '&', sdata, '&', '...', $
               '&', '...', '&', result2[0], '&', 'PL', '&', '...', '\\', $
                format='(a,a,f7.2,a,e11.2,a,e11.2,a,f7.2,a,f7.2,a, a ,a, a ,a, f7.2 ,a, a ,a, a,a)'
    endif 
    if (finite(result1[0]) and finite(result2[0]) and dof1 ne 0 and dof2 ne 0) then begin
        case (abs(red_chi1-1) gt abs(red_chi2-1)) of
        0: begin 
        ;=========================================for shadows====================================================
        cgpolygon, [fitrange_g[0], fitrange_g[1], fitrange_g[1], fitrange_g[0], fitrange_g[0]], $
          [yrange[0], yrange[0], yrange[1], yrange[1], yrange[0]], /data, color = 'dodger blue',linestyle = 1, thick = 3;, /fill
        ;=======================================================================================================
        cgplot, bincenters, gauss1(bincenters, result1), color = 'dodger blue', /overplot, linestyle = 2, thick = 8
        center = String(result1[0], FORMAT='(F0.2)')
        sigma = String(result1[1], FORMAT='(F0.2)')
        text1 = textoidl('\mu = ')
        text2 = textoidl('\sigma = ')
        cgtext, 0.8, 0.6, /data, text1 + center, COLOR = 'blue', charsize = 1
        cgtext, 0.8, 0.2, /DATA, text2 + sigma, COLOR = 'blue', charsize = 1
        printf, lun, name[i], '&', d[i], '&', mass, '&', N0, '&', mudata, '&', sdata, '&', result1[0], $
               '&', result1[1], '&', '...', '&', 'LN', '&', '...', '\\', $
                format='(a,a,f7.2,a,e11.2,a,e11.2,a,f7.2,a,f7.2,a,f7.2,a,f7.2,a, a ,a, a ,a, a,a)' ;  log-normal fitted
        width = [width, result1[1]]
        mach = [mach, ms]
        end
        
        1: begin 
        ;=========================================for shadows====================================================
        cgpolygon, [fitrange_p[0], fitrange_p[1], fitrange_p[1], fitrange_p[0], fitrange_p[0]], $
          [yrange[0], yrange[0], yrange[1], yrange[1], yrange[0]], /data, color = 'magenta',linestyle = 1, thick = 3;, /fill
        ;=======================================================================================================
        cgplot, bincenters, exp(result2[0]*bincenters+result2[1]), color = 'magenta', /overplot, linestyle = 2, thick = 5
        center = String(result2[0], FORMAT='(F0.2)')
        sigma = String(result2[1], FORMAT='(F0.2)')
        text1 = textoidl('\alpha = ')
        text2 = textoidl('c = ')
        cgtext, 0.8, 0.6, /data, text1 + center, COLOR='magenta', charsize = 1
        printf, lun, name[i], '&', d[i], '&', mass, '&', N0, '&', mudata, '&', sdata, '&', '...', $
               '&', '...', '&', result2[0], '&', 'PL', '&', '...', '\\', $
                format='(a,a,f7.2,a,e11.2,a,e11.2,a,f7.2,a,f7.2,a, a ,a, a ,a, f7.2 ,a, a ,a, a,a)' ; power-law fitted
        end
        endcase
    endif 
endfor
free_lun, lun
cgps_close
end