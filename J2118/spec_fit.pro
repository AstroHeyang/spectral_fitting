;+
; NAME:
;   spec_fit
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;     https://github.com/AstroHeyang
;
; PURPOSE:
;   Fit the spectrum using mpfit. 
;
; CALLING SEQUENCE:
;   spec_fit, index_fit, hb_func=hb_func, [o3_func = o3_func, outfile=outfile,$
;     nlfixed=nlfixed, blfixed=blfixed, \sdss]) 
;
; DESCRIPTION:
;   NULL
;
; INPUT:
;   index_fit - fitting region, e.g., 4700-5100 A
;   hb_func - model to fit hbeta region, e.g., 2gauss, 3gauss or 1gauss+1lorentz
;
; KEYWORD PARAMETERS:
;   o3_func - model to fit [OIII], the default is 1 gauss
;   outfile - output file name
;   nlfixed - narrow line fixed to a specific value, e.g., 500 km/s
;   blfixed - broad line fixed to a specific value, e.g., 1200 km/s
;   sdss - sdss spectrum or not 
;
; OUTPUT:
;	A fits file, {wave, flux, error, ra, dec, z, model_x, model_y, chisq, $
;   reduced_chisq, pl, Hb4861_N, Hb4861_B, OIII4959, OIII5007}
;-

function powerlaw, x, p
  if n_elements(p) lt 2 then begin
    print, "powlaw: n_elements(p) lt 2. ERROR!"
    return, !VALUES.F_NAN
  endif
  
  wave = exp(x)
  y = p[0] * (wave/5000.d)^p[1]
  
  return,y
end

function model_2g, x, p
  if n_elements(p) lt 6 then begin
    print, "model_2g: n_elements(p) lt 6. ERROR!"
    return, !VALUES.F_NAN
  endif
  
  flux  = gauss1(x, [p[0], p[1], p[2]])
  flux1 = flux  + gauss1(x, [p[3], p[4], p[5]])
  return, flux1
end

function model_3g, x, p
  if n_elements(p) lt 9 then begin
    print, "model_3g: n_elements(p) lt 9. ERROR!"
    return, !VALUES.F_NAN
  endif
  flux  = gauss1(x, [p[0], p[1], p[2]])
  flux1 = flux  + gauss1(x, [p[3], p[4], p[5]])
  flux2 = flux1  + gauss1(x, [p[6], p[7], p[8]])
  return, flux2
end

function model_1g1l, x, p
  if n_elements(p) lt 6 then begin
    print, "model_1g1l: n_elements(p) lt 6. ERROR!"
    return, !VALUES.F_NAN
  endif
  
  flux  = gauss1(x, [p[0], p[1], p[2]])
  flux1 = flux  + lorentz1(x, [p[3], p[4], p[5]])
  return, flux1
end

function model_expr, indexes, hb_func
  model = 'powerlaw(X, [P[' + strtrim(string(indexes[0]), 2) + ':' + strtrim(string(indexes[1]-1), 2) + ']]) +' ;powerlaw
  
  model += $
    'model_' + hb_func + '(X, [P[' + strtrim(string(indexes[1]), 2) + ':' + strtrim(string(indexes[2]-1), 2) + ']]) +' ;Hbeta
  
  model += $
    'gauss1(X, [P[' + strtrim(string(indexes[2]), 2) + ':' + strtrim(string(indexes[3]-1), 2) + ']]) +' + $ ;OIII4959
    'gauss1(X, [P[' + strtrim(string(indexes[3]), 2) + ':' + strtrim(string(indexes[4]-1), 2) + ']])' ;OIII5007
  	
  return, model
end

PRO spec_fit, index_fit, hb_func = hb_func, o3_func = o3_func, outfile=outfile,$
  nlfixed=nlfixed, blfixed=blfixed, sdss=sdss

  if keyword_set(sdss) then begin
    spec = read_spec(/sdss)
  endif else begin
    spec = read_spec()
  endelse
  wave = spec.wave
  wave_ln = alog(wave)
  flux = spec.flux
  error = spec.flux_error
  z = spec.z 
  n_pixel = n_elements(wave)
  
  if not keyword_set(index_fit) then begin
    index_fit = [4750, 5250]
  endif
  
  if not keyword_set(hb_func) then begin
    hb_func = '2g'
  endif 
  if not keyword_set(o3_func) then begin
    o3_func = '1g'
  endif 
  
  index_labels = ['index_pl', 'index_hb', 'index_oiii4959', 'index_oiii5007']
  n_funcs = n_elements(index_labels) 
  steps = lonarr(n_funcs)
  
  steps[0] = 2 ; powerlaw
  if hb_func eq '3g' then begin
    steps[1] = 9
  endif else begin
    steps[1] = 6
  endelse
  if o3_func eq '2g' then begin
    steps[2] = 6
    steps[3] = 6
  endif else begin
    steps[2] = 3
    steps[3] = 3
  endelse
  
  if not keyword_set(o3_func) then begin
    o3_func = '1g'
  endif
  
  indexes = lonarr(n_funcs + 1)
  for i = 1l, n_funcs do begin
    indexes[i] = indexes[i-1] + steps[i-1]
  endfor
  
  for i = 0l, n_funcs-1 do begin
    void = execute(index_labels[i] + ' = ' + string(indexes[i]))
  endfor
  
  num_P = indexes[n_elements(indexes)-1]
  model = model_expr(indexes, hb_func)
  	
  info_P_str = { $
    value: 0., $
    fixed: 0., $
    limited: [0, 0], $
    limits: [0., 0.],$
    tied: ''}
  info_P = replicate(info_P_str, num_P)
  
  info_P[index_pl].value = 1 
  info_P[index_pl].limited = [0, 0] 
  info_P[index_pl].limits = [0, 0] 
  info_P[index_pl+1].value = 1 
  info_P[index_pl+1].limited = [0, 0] 
  info_P[index_pl+1].limits = [-5, 5]
  
  ; narrow Hbeta
  info_P[index_hb].value = alog(4862.68)
  info_P[index_hb].limited = [1, 1] 
  info_P[index_hb].limits = [alog(4832.), alog(4892.)]
  
  if keyword_set(nlfixed) then begin
    info_P[index_hb+1].value = nlfixed/(2.354*3D5)
    info_P[index_hb+1].fixed = 1  
  endif else begin
  	info_P[index_hb+1].tied =   $
      'P['+ $
      strcompress(string(index_oiii5007+1), /remove_all)+ $
      ']'
  endelse
  info_P[index_hb+2].value = 100./info_P[index_hb].value
  info_P[index_hb+2].limited = [1, 0] 
  info_P[index_hb+2].limits = [0, 0]
  index_bhb1 = index_hb + 3
  
  ; broad Hbeta
  info_P[index_bhb1].value = alog(4872.68)
  info_P[index_bhb1].limited = [1, 1] 
  info_P[index_bhb1].limits = [alog(4832.), alog(4992.)]
  
  if keyword_set(blfixed) then begin
    info_P[index_bhb1+1].value = blfixed/(2.354*3D5)
    info_P[index_bhb1+1].fixed = 1 
  endif else begin 
    info_P[index_bhb1+1].value = 2000./(3D+5*2.354)
    info_P[index_bhb1+1].limited = [1, 1] 
    info_P[index_bhb1+1].limits = [0,3000]/(3D+5*2.354)
  endelse
  info_P[index_bhb1+2].value = 100./info_P[index_hb].value
  info_P[index_bhb1+2].limited = [1, 0] 
  info_P[index_bhb1+2].limits = [0, 0]
  index_bhb2 = index_bhb1 + 3
  
  if hb_func eq '3g' then begin
    info_P[index_bhb2].value = alog(4862.68)
    info_P[index_bhb2].limited = [1, 1] 
    info_P[index_bhb2].limits = [alog(4832.), alog(4922.)]
    info_P[index_bhb2+1].value = 2000./(3D+5*2.354)
    info_P[index_bhb2+1].limited = [1, 1] 
    info_P[index_bhb2+1].limits = [800, 7000]/(3D+5*2.354)
    info_P[index_bhb2+2].value = 100./info_P[index_hb].value
    info_P[index_bhb2+2].limited = [1, 0] 
    info_P[index_bhb2+2].limits = [0, 0]
  endif
  
  info_P[index_oiii5007].value		= alog(5008.24)
  info_P[index_oiii5007].limited		= [1, 1] 
  info_P[index_oiii5007].limits		= [alog(4980.), alog(5040.)]
  info_P[index_oiii5007+1].value		= 300./(3D+5*2.354)
  info_P[index_oiii5007+1].limited	= [1, 1] 
  info_P[index_oiii5007+1].limits		= [0, 800]/(3D+5*2.354)
  info_P[index_oiii5007+2].value		= 100./info_P[index_oiii5007].value
  info_P[index_oiii5007+2].limited	= [1, 0] 
  info_P[index_oiii5007+2].limits		= [0, 0]
  
  info_P[index_oiii4959].tied   =   $
    'P['+ $
    strcompress(string(index_oiii5007), /remove_all)+ $
    '] - alog(5008.24/4960.30)'
  info_P[index_oiii4959+1].tied   =   $
    'P['+ $
    strcompress(string(index_oiii5007+1), /remove_all)+ $
    ']'
  info_P[index_oiii4959+2].tied   =   $
    '1./2.98 * P['+ $
    strcompress(string(index_oiii5007+2), /remove_all)+ $
    ']'
  
  if o3_func eq '2g' then begin
  
  	info_P[index_oiii5007+3].value		= 5008.24
  	info_P[index_oiii5007+3].limited	= [1, 1] 
  	info_P[index_oiii5007+3].limits		= [4978, 5038]
  	info_P[index_oiii5007+4].value		= 2000./(3D+5*2.354)
  	info_P[index_oiii5007+4].limited	= [1, 0] 
  	info_P[index_oiii5007+4].limits		= [0, 0]
  	info_P[index_oiii5007+5].value		= 100./info_P[index_oiii5007].value
  	info_P[index_oiii5007+5].limited	= [1, 0] 
  	info_P[index_oiii5007+5].limits		= [0, 0]
  
  	info_P[index_oiii4959+3].tied   =   $
      'P['+ $
      strcompress(string(index_oiii5007+3), /remove_all)+ $
      '] - alog(6564.61/4862.68)'
  
  	info_P[index_oiii4959+4].tied   =   $
      'P['+ $
      strcompress(string(index_oiii5007+4), /remove_all)+ $
      ']'

  	info_P[index_oiii4959+5].tied   =   $
      '1./3. * P['+ $
      strcompress(string(index_oiii5007+5), /remove_all)+ $
      ']'
  
  endif
  
  ind_fit = where(wave ge index_fit[0] and wave le index_fit[1])
  fit_x = wave_ln[ind_fit]
  fit_y = flux[ind_fit]
  fit_error = error[ind_fit]
  
  time0   = systime(1)
  params  = mpfitexpr( $
    model, fit_x, fit_y, fit_err, $
    Parinfo=info_P, $
    yfit=yfit, perror=perror, $
    bestnorm=bestnorm, dof=dof, $
    status=status,maxiter=500,gtol=1.d-10,$
    ftol=1.d-10,/quiet)
  time1   = systime(1)
  print, ''
  print, 'Fitting Time Consumed: ', time1-time0, 'Seconds', $
    format='(A25, D20.3, A8)'
  print, 'status:', format='(A, $)'
  print, status
  print, 'Chisq:', format='(A, $)'
  print, bestnorm, dof, bestnorm/dof
  print, 'Chisq calculated manually:'
  chisq = (yfit - fit_y)^2/(fit_error^2)
  n_dofs = n_elements(fit_y) -  (2+6+3-1)
  print, total(chisq), n_dofs, total(chisq)/n_dofs
  
  ;output
  out_str = {					$
    wave: wave,				$
    flux: flux,				$
    error: error,			$
    ra: spec.ra,			$
    dec: spec.dec,			$
    z: spec.z,				$
    model_x: exp(fit_x),	$
    model_y: yfit,			$
    chisq: dblarr(n_elements(fit_y)),	$
    reduced_chisq:0.d,		$
    pl: dblarr(4),			$
    Hb4861_N: dblarr(6),	$
    Hb4861_B: dblarr(12),	$
    OIII4959: dblarr(12), 	$
    OIII5007: dblarr(12)  	$
  	}
  
  output = replicate(out_str, 1)
  output.chisq = chisq
  output.reduced_chisq = total(chisq)/n_dofs
  
  output.pl[0] = params[index_pl]
  output.pl[1] = perror[index_pl]
  output.pl[2] = params[index_pl+1]
  output.pl[3] = perror[index_pl+1]
  
  for i=0,2 do begin
    output.Hb4861_N[2*i] =   params[index_hb + i]
    output.Hb4861_N[2*i+1] = perror[index_hb + i]
    output.Hb4861_B[2*i] =   params[index_hb + 3 +i]
    output.Hb4861_B[2*i+1] = perror[index_hb + 3 + i]
    
    if hb_func eq '3g' then begin
      output.Hb4861_B[2*(3+i)] =   params[index_hb + 6 +i]
      output.Hb4861_B[2*(3+i)+1] = perror[index_hb + 6 + i]
    endif
    
    output.OIII4959[2*i]	= params[index_oiii4959 + i]
    output.OIII4959[2*i+1]	= perror[index_oiii4959 + i]
    output.OIII5007[2*i]	= params[index_oiii5007 + i]
    output.OIII5007[2*i+1]	= perror[index_oiii5007 + i]
    
    if o3_func eq '2g' then begin
      output.OIII4959[2*(3+i)]	= params[index_oiii4959 + i + 6]
      output.OIII4959[2*(3+i)+1]= perror[index_oiii4959 + i + 6]
      output.OIII5007[2*(3+i)]	= params[index_oiii5007 + i + 6]
      output.OIII5007[2*(3+i)+1]= perror[index_oiii5007 + i + 6]
    endif
  endfor
  
  if not keyword_set(outfile) then begin
    outfile = 'line_J2118.fits'
  endif
  
  mwrfits, output, outfile, /create
  print, '0O0 Finished! 0O0'
END
