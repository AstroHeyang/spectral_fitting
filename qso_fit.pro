;+
; NAME:
;   nlFeIIfunc
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;       https://github.com/AstroHeyang
;
; PURPOSE:
;   fit the optical spectrum of a qso
;
; CALLING SEQUENCE:
;   qso_fit, 
;
; DESCRIPTION:
;   The Fe II multiplets are obtained from the well-known NLS1 Izw~1 (Veron-City et al. 2004)
;   Also see Dong et al. (2011) for details.
;
; INPUT:
;   x - an array of rest-wavelength in log-scale
;
;   p - an array including 3 parameters of Fe II model;
;       [log_w0, sigma, flux] for gauss1 function;
;
; KEYWORD PARAMETERS:
;   forbidden - only forbidden Fe II lines are used
;
;   permitted - only permitted Fe II lines are used
;
;   extended - forbidden and permitted Fe II lines and other lines
;              including Cr, Ni, Si, Ti, Ca, are used

; OUTPUT:
;   return the array of flux for narrow-line Fe II
;
; VERSIONS:
;   v1 - xbdong 2005
;-

function read_ascii_file, spec_file
	readcol, spec_file, wave, flux, flux_err, format = 'D, D, D'
	spec = {				$
		wave: wave,			$
		flux: flux,			$
		flux_err: flux_err	$
		}
	return spec
end

function get_func_index, par
; get index for each component in the  mpfit function 
	func_index = {
		broken_powerlaw:	[0, 3],	$
		nl_FeII_func:		[0, 3],	$ 
		bl_FeII_func:		[0, 3],	$
        nl_Hgamma_4340:		[0, 3],	$
        bl_Hgamma_4340:		[0, 3],	$
        OIII_4363:			[0, 3],	$
		nl_HeII_4686:		[0, 3],	$     
        bl_HeII_4686:		[0, 3],	$
        nl_Hbeta_4861:		[0, 3],	$
        bl_Hbeta_4861:		[0, 3],	$
        OIII_4959:			[0, 3],	$	
        OIII_5007:			[0, 3],	$	
        NI_5199:			[0, 3],	$	
        nl_HeI_5876:		[0, 3],	$	
        bl_HeI_5876:		[0, 3],	$	
		OI_6300:			[0, 3],	$	
		SIII_6313:			[0, 3],	$	
		OI_6364:			[0, 3],	$	
		NII_6548:			[0, 3],	$	
		nl_Halpha_6563:		[0, 3],	$
		bl_Halpha_6563:		[0, 3],	$
		NII_6583:			[0, 3],	$	
		HeI_6678:			[0, 3],	$	
        SII_6716:			[0, 3],	$	
        SII_6731:			[0, 3],	
		}

	func_index_tags = tag_names(func_index)
	par_tags = tag_names(par)
	for i = 2, n_elements(func_index_tags)-1 do begin
		index_this_tag = where(par_tags eq func_index_tags[i])
		;; current_index = last_index + nGauss * 3, nGauss is obtained from 
		;; the par file, else nGauss = 1 (default)
		func_index.(i)[0] = func_index.(i-1)[0] + func_index.(i-1)[1]
		if index_this_tag[0] ne -1 then begin
			func_index.(i)[1] = long(par.(index_this_tag)) * 3
		endif 
	endfor
	return, func_index
end

function get_mpfit_func, line_infos
;; generate function for mpfit
	nComponents = n_elements(line_infos)
	func = 'broken_powerlaw(X, [P[0:2]]) +'			 $
			+ 'nlFeIIfunc(X, [P[3:5]], /extended) +' $
			+ 'blFeIIfunc(X, [P[6:8]], /extended) +' $
	for i = 3, nComponents-1 do begin
		if i lt nComponents-1 then begin
			func += 'ngauss(X, [P[' + strtrim(string(func_index.(i)[0]), 2) + ':' + $
					strtrim(string(func_index.(i)[0] + func_index.(i)[1]-1), 2) + ']]) + '
		endif else begin
			func += 'ngauss(X, [P[' + strtrim(string(func_index.(i)[0]), 2) + ':' + $
					strtrim(string(func_index.(i)[0] + func_index.(i)[1]-1), 2) + ']]) '
		endelse
	endfor

	return, func
end

function get_parameter_info, info_P, func_index, fixed2OIII=fixed2OIII, $
	fixed2SII=fixed2SII, fixed2NII=fixed2NII
;; get informations for parameters, such as initial values, constraints
;; these parameters are set for mpfit
;; here we use the method of exhaustion

	line_infos = get_line_info()
	n_lines = n_elements(line_infos)

	;; set initial values
	line_tags = tag_names(func_index)
	for i = 0, n_elements(line_tags)-1 do begin
		this_index = func_index.(i)[0]
		n_parameter = func_index.(i)[1]
		line_index = where(strupcase(line_infos.name) eq line_tags[i])
		this_values = line_infos[line_index].values
		;; TODO: change the initial values for nGauss model
		while n_parameter ge 3 do begin
			info_P[this_index:this_index+2].value = this_values
			n_parameter -= 3
			this_index += 3
		endwhile 
			
		while this_index lt next_index:
			info_P[this_index:this_index+2].value = this_values
			this_index += 3
		endwhile
	endfor

	;; set tied, limits and limited
	for i = 0, n_elements(line_tags)-1 do begin
		this_index = func_index.(i)[0]
		

		line_index = where(strupcase(line_infos.name) eq line_tags[i])
		this_tied = line_infos[line_index].tied
		for j = 0,n_elements(this_tied) do begin

			if this_tied[j] ne -1 then begin
				line_index = line_infos[this_tied[j]]
				tied_index = line_infos.name[]
			info_P[this_index].tied = $
				'P[' + strcompress
	
			;; if no tied, then limited and limits are set
	

	index_nl_model = line_tags.(where(line_tags eq nl_model))
	index_bl_model = line_tags.bl_Halpha_6563 ; default: broad lines fixed to broad Ha
	for i = 0, n_elements(line_tags)-2 do begin
		this_index = func_index.(i)
		next_index = func_index.(i+1)
		case line_tags[i] of 
			strupcase('broken_powerlaw'): $
				info_P[this_index:this_index+2].limited = [[1, 0], [1,1], [1,1]]
				info_P[this_index:this_index+2].limits = [[0., 0.], [-5., 5.], [-5., 5.]]
			strupcase('OIII_5007'): $
				
			else: begin
				if initial_values.(i)[4] eq 0 then begin
					info_P[this_index].tied = $
						'P[' + strcompress(string(nl_centroid_index), /remove_all) + $
						'] - ' + strcompress(string(initial_values.(index_nl_model)[0]), $
						/remove_all)) + '+' + strcompress(string(initial_values.(i)[0]), $
						/remove_all))
					info_P[this_index+1].tied = $
						'P[' + strcompress(string(index_nl_model), /remove_all) + $
						'] - ' + strcompress(string(initial_values.(index_nl_model)[0]), $
						/remove_all)) + '+' + strcompress(string(initial_values.(i)[0]), $
						/remove_all))
				else begin


	
        info_P[ind_hb4861_b3].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b3), /remove_all)+ $
                                        '] - alog(6564.61/4862.68)'
        info_P[ind_hb4861_b3+1].tied =   'P['+ $
                                        strcompress(string(ind_ha6563_b3+1), /remove_all)+ $
                                        ']'
        info_P[ind_hb4861_b3+2].limited= [1, 0]
        info_P[ind_hb4861_b3+2].limits=  [0., 0.]

		if lines[i] eq strupcase('broken_powerlaw') then begin
			ind_pl = func_index.(i)
			info_P[ind_pl].value    = 1.
			info_P[ind_pl+1].value	= -1.7
			info_P[ind_pl+2].value  = -1.7
			info_P[ind_pl].limited  =   [1, 0]
			info_P[ind_pl].limits   =   [0., 0.]
			info_P[ind_pl+1].limited=   [1, 1]
			info_P[ind_pl+1].limits =   [-5., 5.]
			info_P[ind_pl+2].limited=   [1, 1]
			info_P[ind_pl+2].limits =   [-5., 5.]
		endif

    ; narrow FeII lines
    ; width tied to [OIII]5007, flux constrainted
    info_P[ind_nl_FeII+1].tied   =   'P['+ $
                                    strcompress(string(ind_oiii5007+1), /remove_all)+ $
                                    ']'
    info_P[ind_nl_FeII+2].limited=   [1, 0]
    info_P[ind_nl_FeII+2].limits =   [0., 0.]

PRO qso_fit, wave, flux, z, flux_err=flux_err, ebv=ebv, ra=ra, dec=dec  

	On_error, 2
	par = read_par('default.par')
	FeII_file = par.FeII_file
	output_dir = par.output_dir 
    if not file_test(out_dir, /dir) then begin
        file_mkdir, out_dir
    endif

	;; redshift correction
	log_wave = alog(wave)
	rest_log_wave = log_wave - alog(1.0 + z)
	rest_wave = exp(rest_log_wave)

	;; Galactic extinction correction
	if not keyword_set(ebv) then begin
		glactc, spec_ra, spec_dec, 2000, gal_l, gal_b, 1, /degree
		ebv = dust_getval(gal_l, gal_b, ipath='./SFD/maps/',/interp)
	endif 
    fm_unred, rest_wave, flux, ebv, unred_flux
    fm_unred, rest_wave, flux_err, ebv, unred_flux_err

	;; get spectra with good pixels
	good_pixels = where(flux_err gt 0)
	if good_pixels[0] ne -1 then begin
		rest_log_wave = rest_log_wave[good_pixels]
		unred_flux = unred_flux[good_pixels]
		unred_flux_err = unred_flux_err[good_pixels] 
	else begin
		print('ERROR! No good pixels!')
		stop
	endelse

;===============================================================================
	;; parameters
	func_index = get_func_index(par)
	model_func = get_mpfit_func(func_index)

	num_P = func_index.num_P + 1
    info_P_str = {  value:      0., $
                    limited:    [0, 0], $
                    limits:     [0., 0.],$
                    tied:       ''}
    info_P = replicate(info_P_str, num_P)

	;;;;;;;;;;;;;;;;;;;;;;;; broad line centroid wavelength and width, you could revised these to obtain better results;;;;;;;;;

    ;ha6563b_wave = [6564.61,6584.61,6594.61,6634.61]
    hb4861b_wave = [4792.68,4802.68,4862.68,4892.68]
    ha6563b_wave = [6544.61,6564.61,6584.61,6594.61]
    ;hb4861b_wave = [4862.68,4862.68,4862.68,4862.68]
    ha6563b_width = [3000,10000,5000,3000]
    hb4861b_width = [3000,10000,5000,3000]
	   if bl eq '3g' then begin
        info_P[ind_hc4340_b3].value  = alog(4341.70)
        info_P[ind_hb4861_b3].value  = alog(hb4861b_wave[2])
        info_P[ind_ha6563_b3].value  = alog(ha6563b_wave[2])
        info_P[ind_hc4340_b3+1].value    = 500./(3D+5*2.354)
        info_P[ind_hb4861_b3+1].value   = ha6563b_width[2]/(3D+5*2.354)
        info_P[ind_ha6563_b3+1].value   = hb4861b_width[2]/(3D+5*2.354)
        info_P[ind_hc4340_b3+2].value    = 100./info_P[ind_hc4340_b3].value
        info_P[ind_hb4861_b3+2].value    = 500./info_P[ind_hb4861_b3].value
        info_P[ind_ha6563_b3+2].value    = 500./info_P[ind_ha6563_b3].value
        info_P[ind_ha6563_b3].limited    =   [1, 1]
        ;info_P[ind_ha6563_b3].limits     =   [alog(6533), alog(6593)]
        info_P[ind_ha6563_b3].limits     =   [alog(6403), alog(6653)]
    ;   info_P[ind_ha6563_b3].limits     =   [alog(6303), alog(6753)]
        info_P[ind_ha6563_b3+1].limited  =   [1, 1]
        info_P[ind_ha6563_b3+1].limits   =   [400, 2d4]/(3d5*2.354)
        info_P[ind_ha6563_b3+2].limited  =   [1, 0]
        info_P[ind_ha6563_b3+2].limits   =   [0., 0.]

        info_P[ind_hb4861_b3].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b3), /remove_all)+ $
                                        '] - alog(6564.61/4862.68)'
        info_P[ind_hb4861_b3+1].tied =   'P['+ $
                                        strcompress(string(ind_ha6563_b3+1), /remove_all)+ $
                                        ']'
        info_P[ind_hb4861_b3+2].limited= [1, 0]
        info_P[ind_hb4861_b3+2].limits=  [0., 0.]

        info_P[ind_hc4340_b3].tied       =   'P['+ $
                                        strcompress(string(ind_ha6563_b3), /remove_all)+ $
                                        '] - alog(6564.61/4341.68)'
        info_P[ind_hc4340_b3+1].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b3+1), /remove_all)+ $
                                        ']'
        info_P[ind_hc4340_b3+2].limited  =   [1, 0]
        info_P[ind_hc4340_b3+2].limits   =   [0., 0.]
    endif

	    if bl eq '4g' then begin
        info_P[ind_hc4340_b3].value  = alog(4341.70)
        info_P[ind_hc4340_b4].value  = alog(4341.70)
        info_P[ind_hb4861_b3].value  = alog(hb4861b_wave[2])
        info_P[ind_hb4861_b4].value  = alog(hb4861b_wave[3])
        info_P[ind_ha6563_b3].value  = alog(ha6563b_wave[2])
        info_P[ind_ha6563_b4].value  = alog(ha6563b_wave[3])
        info_P[ind_hc4340_b3+1].value   = 500./(3D+5*2.354)
        info_P[ind_hc4340_b4+1].value   = 500./(3D+5*2.354)
        info_P[ind_hb4861_b3+1].value   = hb4861b_width[2]/(3D+5*2.354)
        info_P[ind_hb4861_b4+1].value   = hb4861b_width[3]/(3D+5*2.354)
        info_P[ind_ha6563_b3+1].value   = ha6563b_width[2]/(3D+5*2.354)
        info_P[ind_ha6563_b4+1].value   = ha6563b_width[3]/(3D+5*2.354)
        info_P[ind_hc4340_b3+2].value    = 200./info_P[ind_hc4340_b3].value
        info_P[ind_hc4340_b4+2].value    = 200./info_P[ind_hc4340_b4].value
        info_P[ind_hb4861_b3+2].value    = 500./info_P[ind_hb4861_b3].value
        info_P[ind_hb4861_b4+2].value    = 500./info_P[ind_hb4861_b4].value
        info_P[ind_ha6563_b3+2].value    = 500./info_P[ind_ha6563_b3].value
        info_P[ind_ha6563_b4+2].value    = 500./info_P[ind_ha6563_b4].value
        info_P[ind_ha6563_b3].limited    =   [1, 1]
    ;   info_P[ind_ha6563_b3].limits     =   [alog(6433), alog(6593)]
        info_P[ind_ha6563_b3].limits     =   [alog(6303), alog(6753)]
        info_P[ind_ha6563_b3+1].limited  =   [1, 1]
        info_P[ind_ha6563_b3+1].limits   =   [400, 1d4]/(3d5*2.354)
        info_P[ind_ha6563_b3+2].limited  =   [1, 0]
        info_P[ind_ha6563_b3+2].limits   =   [0., 0.]

        info_P[ind_ha6563_b4].limited    =   [1, 1]
        info_P[ind_ha6563_b4].limits     =   [alog(6303), alog(6753)]
        info_P[ind_ha6563_b4+1].limited  =   [1, 1]
        info_P[ind_ha6563_b4+1].limits   =   [300, 1d4]/(3d5*2.354)
        info_P[ind_ha6563_b4+2].limited  =   [1, 0]
        info_P[ind_ha6563_b4+2].limits   =   [0., 0.]
		info_P[ind_hb4861_b3].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b3), /remove_all)+ $
                                        '] - alog(6564.61/4862.68)'
        info_P[ind_hb4861_b3+1].tied =   'P['+ $
                                        strcompress(string(ind_ha6563_b3+1), /remove_all)+ $
                                        ']'
        info_P[ind_hb4861_b3+2].limited= [1, 0]
        info_P[ind_hb4861_b3+2].limits=  [0., 0.]

        info_P[ind_hb4861_b4].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b4), /remove_all)+ $
                                        '] - alog(6564.61/4862.68)'
        info_P[ind_hb4861_b4+1].tied =   'P['+ $
                                        strcompress(string(ind_ha6563_b4+1), /remove_all)+ $
                                        ']'
        info_P[ind_hb4861_b4+2].limited= [1, 0]
        info_P[ind_hb4861_b4+2].limits=  [0., 0.]



        info_P[ind_hc4340_b3].tied       =   'P['+ $
                                        strcompress(string(ind_ha6563_b3), /remove_all)+ $
                                        '] - alog(6564.61/4341.68)'
        info_P[ind_hc4340_b3+1].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b3+1), /remove_all)+ $
                                        ']'
        info_P[ind_hc4340_b3+2].limited  =   [1, 0]
        info_P[ind_hc4340_b3+2].limits   =   [0., 0.]
        info_P[ind_hc4340_b4].tied       =   'P['+ $
                                        strcompress(string(ind_ha6563_b4), /remove_all)+ $
                                        '] - alog(6564.61/4341.68)'
        info_P[ind_hc4340_b4+1].tied   =   'P['+ $
                                        strcompress(string(ind_ha6563_b4+1), /remove_all)+ $
                                        ']'
        info_P[ind_hc4340_b4+2].limited  =   [1, 0]
        info_P[ind_hc4340_b4+2].limits   =   [0., 0.]
    endif

	;fit parameters
    info_P[ind_hc4340_n].value  = alog(4341.68)
    info_P[ind_hc4340_b1].value  = alog(4341.68)
    info_P[ind_hc4340_b2].value  = alog(4341.68)
    info_P[ind_oiii4363].value  = alog(4364.44)
    info_P[ind_nl_FeII].value    = alog(5019.83)
    info_P[ind_bl_FeII].value    = alog(5019.85)
    info_P[ind_heii4686_n].value= alog(4686.)
    info_P[ind_heii4686_b].value= alog(4687.)
    info_P[ind_hb4861_n].value  = alog(4861.)
    info_P[ind_hb4861_b1].value  = alog(4862.)
    info_P[ind_hb4861_b2].value  = alog(4862.)
    info_P[ind_oiii4959].value  = alog(4959.)
    info_P[ind_oiii5007].value  = alog(5007.)
    ;if ( max(good_wave) gt 7000. ) then begin
    info_P[ind_hei5876_n].value = alog(5876.)
    info_P[ind_hei5876_b].value = alog(5877.)
    info_P[ind_oi6300].value    = alog(6300.)
    info_P[ind_siii6313].value  = alog(6313.)
    info_P[ind_oi6364].value    = alog(6364.)
    info_P[ind_nii6548].value   = alog(6548.)
    info_P[ind_nii6583].value   = alog(6583.)
    info_P[ind_ha6563_n].value  = alog(6563.)
    info_P[ind_ha6563_b1].value  = alog(6564.)
    info_P[ind_ha6563_b2].value  = alog(6564.)
    info_P[ind_sii6716].value   = alog(6716.)
    info_P[ind_sii6731].value   = alog(6731.)
    ;endif

    ; lines width

    info_P[ind_hc4340_n+1].value    = 500./(3D+5*2.354)
    info_P[ind_hc4340_b1+1].value    = 2000./(3D+5*2.354)
    info_P[ind_hc4340_b2+1].value    = 2000./(3D+5*2.354)
    info_P[ind_oiii4363+1].value    = 500./(3D+5*2.354)
    info_P[ind_nl_FeII+1].value      = 500./(3D+5*2.354)
    info_P[ind_bl_FeII+1].value      = 3000./(3D+5*2.354)
    info_P[ind_heii4686_n+1].value  = 500./(3D+5*2.354)
    info_P[ind_heii4686_b+1].value  = 3000./(3D+5*2.354)
	info_P[ind_oiii4959+1].value    = 50./(3D+5*2.354)
    info_P[ind_oiii5007+1].value    = 50./(3D+5*2.354)
    ;if ( max(good_wave) gt 7000. ) then begin
    info_P[ind_hei5876_n+1].value   = 500./(3D+5*2.354)
    info_P[ind_hei5876_b+1].value   = 2000./(3D+5*2.354)
    info_P[ind_oi6300+1].value      = 500./(3D+5*2.354)
    info_P[ind_siii6313+1].value    = 500./(3D+5*2.354)
    info_P[ind_oi6364+1].value      = 500./(3D+5*2.354)
    info_P[ind_nii6548+1].value     = 300./(3D+5*2.354)
    info_P[ind_nii6583+1].value     = 300./(3D+5*2.354)
    info_P[ind_ha6563_n+1].value    = 300./(3D+5*2.354)
    info_P[ind_ha6563_b1+1].value   = 2000./(3D+5*2.354)
    info_P[ind_ha6563_b2+1].value   = 2000./(3D+5*2.354)
    info_P[ind_sii6716+1].value     = 300./(3D+5*2.354)
    info_P[ind_sii6731+1].value     = 300./(3D+5*2.354)
    ;endif
	info_P[ind_hc4340_n+2].value    = 50./info_P[ind_hc4340_n].value
    info_P[ind_hc4340_b1+2].value    = 50./info_P[ind_hc4340_b1].value
    info_P[ind_hc4340_b2+2].value    = 50./info_P[ind_hc4340_b2].value
    info_P[ind_oiii4363+2].value    = 10./info_P[ind_oiii4363].value
    info_P[ind_nl_FeII+2].value      = 10./info_P[ind_nl_FeII].value
    info_P[ind_bl_FeII+2].value      = 10./info_P[ind_bl_FeII].value
    info_P[ind_heii4686_n+2].value  = 10./info_P[ind_heii4686_n].value
    info_P[ind_heii4686_b+2].value  = 10./info_P[ind_heii4686_b].value
    info_P[ind_hb4861_n+2].value    = 500./info_P[ind_hb4861_n].value
    info_P[ind_hb4861_b1+2].value    = 1000./info_P[ind_hb4861_b1].value
    info_P[ind_hb4861_b2+2].value    = 1000./info_P[ind_hb4861_b2].value
    info_P[ind_oiii4959+2].value    = 100./info_P[ind_oiii4959].value
    info_P[ind_oiii5007+2].value    = 300./info_P[ind_oiii5007].value
    ;if ( max(good_wave) gt 7000. ) then begin
    info_P[ind_hei5876_n+2].value   = 50./info_P[ind_hei5876_n].value
    info_P[ind_hei5876_b+2].value   = 50./info_P[ind_hei5876_b].value
    info_P[ind_oi6300+2].value      = 50./info_P[ind_oi6300].value
    info_P[ind_siii6313+2].value    = 50./info_P[ind_siii6313].value
    info_P[ind_oi6364+2].value      = 50./info_P[ind_oi6364].value
    info_P[ind_nii6548+2].value     = 200./info_P[ind_nii6548].value
    info_P[ind_nii6583+2].value     = 200./info_P[ind_nii6583].value
    info_P[ind_ha6563_n+2].value    = 2000./info_P[ind_ha6563_n].value
    info_P[ind_ha6563_b1+2].value    = 1000./info_P[ind_ha6563_b1].value
    info_P[ind_ha6563_b2+2].value    = 1000./info_P[ind_ha6563_b2].value
    info_P[ind_sii6716+2].value     = 400./info_P[ind_sii6716].value
    info_P[ind_sii6731+2].value     = 400./info_P[ind_sii6731].value

	; powerlaw continuum
    ;info_P[ind_pl+3].value          = -1.7

    if keyword_set(o32g) then begin
       ; info_P[ind_oiii4959+3].value        = alog(4959.)
        ;info_P[ind_oiii5007+3].value        = alog(5007.)
        info_P[ind_oiii4959+3].value        = alog(4955.)
        info_P[ind_oiii5007+3].value        = alog(5002.)
        info_P[ind_oiii4959+1+3].value      = 100./(3D+5*2.354)
        info_P[ind_oiii5007+1+3].value      = 100./(3D+5*2.354)
        info_P[ind_oiii4959+2+3].value      = 300./info_P[ind_oiii4959].value
        info_P[ind_oiii5007+2+3].value      = 900./info_P[ind_oiii5007].value
    endif

	; [OIII]5007
    ; line center, width and flux constrainted
    info_P[ind_oiii5007].limited    =   [1, 1]
    info_P[ind_oiii5007].limits     =   [alog(4980.), alog(5040.)]
    info_P[ind_oiii5007+1].limited  =   [1, 1]
    ;info_P[ind_oiii5007+1].limits   =   [0.,2000.]/(3d5*2.354)
    info_P[ind_oiii5007+1].limits   =   [0.,800.]/(3d5*2.354)
    info_P[ind_oiii5007+2].limited  =   [1, 0]
    info_P[ind_oiii5007+2].limits   =   [0., 0.]

    if keyword_set(o32g) then begin

        info_P[ind_oiii5007+3].limited  =   [1, 1]
        info_P[ind_oiii5007+3].limits   =   [alog(4970.), alog(5040.)]
        info_P[ind_oiii5007+4].limited  =   [1, 1]
        info_P[ind_oiii5007+4].limits   =   [0.,2000.]/(3d5*2.354)
        info_P[ind_oiii5007+5].limited  =   [1, 0]
        info_P[ind_oiii5007+5].limits   =   [0., 0.]

		endif
    ; [OIII]4959
    ; redshift, width and flux tied to [OIII]5007
    info_P[ind_oiii4959].tied   =   'P['+ $
                                    strcompress(string(ind_oiii5007), /remove_all)+ $
                                    '] - alog(5008.24/4960.30)'
    info_P[ind_oiii4959+1].tied =   'P['+ $
                                    strcompress(string(ind_oiii5007+1), /remove_all)+ $
                                    ']'
    info_P[ind_oiii4959+2].tied =   '1./3. * P['+ $
                                    strcompress(string(ind_oiii5007+2), /remove_all)+ $
                                    ']'
    if keyword_set(o32g) then begin
        info_P[ind_oiii4959+3].tied =   'P['+ $
                                        strcompress(string(ind_oiii5007+3), /remove_all)+ $
                                        '] - alog(5008.24/4960.30)'
        info_P[ind_oiii4959+4].tied =   'P['+ $
                                        strcompress(string(ind_oiii5007+4), /remove_all)+ $
                                        ']'
        info_P[ind_oiii4959+5].tied =   '1./2.98 * P['+ $
                                        strcompress(string(ind_oiii5007+5), /remove_all)+ $
                                        ']'
    endif
	; [OIII]4363
    ; redshift and width tied to [OIII]5007, flux constrainted
    info_P[ind_oiii4363].tied =     'P['+ $
                                    strcompress(string(ind_oiii5007), /remove_all)+ $
                                    '] - alog(5008.24/4364.44)'
    info_P[ind_oiii4363+1].tied=    'P['+ $
                                    strcompress(string(ind_oiii5007+1), /remove_all)+ $
                                    ']'
    info_P[ind_oiii4363+2].limited= [1, 0]
    info_P[ind_oiii4363+2].limits=  [0., 0.]


    ; narrow Ha 6563
    ; line center, width tied to [OIII]5007, [SII] 6731 or set to be free
    info_P[ind_ha6563_n].limited    =   [1, 1]
    info_P[ind_ha6563_n].limits     =   [alog(6533), alog(6593)]


    info_P[ind_ha6563_n+2].limited  =   [1, 0]
    info_P[ind_ha6563_n+2].limits   =   [0., 0.]

	ind_nlinemodel = ind_ha6563_n
    if keyword_set(fixed2oiii) then begin
        ;info_P[ind_ha6563_n+1].tied=    'P['+ $
        ;                                strcompress(string(ind_oiii5007+1), /remove_all)+ $
        ;                                ']'
        info_P[ind_ha6563_n+1].tied=    'P['+ $
                                        strcompress(string(ind_hb4861_n+1), /remove_all)+ $
                                        ']'
        ind_nlinemodel = ind_oiii5007
    endif else if keyword_set(fixed2sii) then begin
        info_P[ind_ha6563_n+1].tied =   'P['+ $
                                        strcompress(string(ind_sii6731+1), /remove_all)+ $
                                        ']'
        ind_nlinemodel = ind_sii6731
    endif else begin
        info_P[ind_ha6563_n+1].limited  =   [1, 1]
        ;info_P[ind_ha6563_n+1].limits   =   [0., 1000.]/(3d5*2.354)
        info_P[ind_ha6563_n+1].limits   =   [0., 450.]/(3d5*2.354)
    endelse

	; narrow Hbeta 4861
    ; redshift tied to narrow Ha 6563, width tied to [OIII]5007 or NHa,flux constrainted
    info_P[ind_hb4861_n].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/4862.68)'
    if keyword_set(nhbfixedoiii) or keyword_set(fixed2oiii) then begin
        info_P[ind_hb4861_n+1].tied =   'P['+ $
                                        strcompress(string(ind_oiii5007+1), /remove_all)+ $
                                        ']'
    endif else begin
        info_P[ind_hb4861_n+1].tied =   'P['+ $
                                        strcompress(string(ind_ha6563_n+1), /remove_all)+ $
                                        ']'
    endelse

	info_P[ind_hb4861_n+2].limited= [1, 0]
    info_P[ind_hb4861_n+2].limits=  [0., 0.]

    ; narrow H_gamma 4340
    ; redshift and width tied to narrow Ha 6563, flux constrainted
    info_P[ind_hc4340_n].tied =     'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/4341.68)'
    info_P[ind_hc4340_n+1].tied=    'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_hc4340_n+2].limited= [1, 0]
    info_P[ind_hc4340_n+2].limits=  [0., 0.]

    ; narrow HeII 4686
    ; redshift and width tied to narrow Ha 6563, flux constrainted
    info_P[ind_heii4686_n].tied =   'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/4687.02)'
    info_P[ind_heii4686_n+1].tied=  'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_heii4686_n+2].limited= [1, 0]
    info_P[ind_heii4686_n+2].limits=  [0., 0.]

	; narrow HeI 5876
    ; redshift and width tied to narrow Ha 6563, flux constrainted
    info_P[ind_hei5876_n].tied =   'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/5877.29)'
    info_P[ind_hei5876_n+1].tied=  'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_hei5876_n+2].limited= [1, 0]
    info_P[ind_hei5876_n+2].limits=  [0., 0.]
    ; narrow OI 6300
    ; redshift tied to narrow Ha 6563, width tied to narrow-line model,flux constrainted
    info_P[ind_oi6300].tied =       'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/6302.05)'
    info_P[ind_oi6300+1].tied=      'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_oi6300+2].limited= [1, 0]
    info_P[ind_oi6300+2].limits=  [0., 0.]

    ; narrow SIII 6313
    ; redshift tied to narrow Ha 6563, width tied to narrow-line model,flux constrainted
    info_P[ind_siii6313].tied =     'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/6313.8)'
    info_P[ind_siii6313+1].tied=    'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_siii6313+2].limited= [1, 0]
    info_P[ind_siii6313+2].limits=  [0., 0.]

	; narrow OI 6364
    ; redshift tied to narrow Ha 6563, width tied to narrow-line model,flux constrainted
    info_P[ind_oi6364].tied =       'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/6365.54)'
    info_P[ind_oi6364+1].tied=      'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_oi6364+2].limited= [1, 0]
    info_P[ind_oi6364+2].limits=  [0., 0.]

    ; [NII] 6548
    ; redshift and width tied to narrow Ha 6563, flux constrainted
    info_P[ind_nii6548].tied    =   'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/6549.85)'
    info_P[ind_nii6548+1].tied  =   'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_nii6548+2].limited=  [1, 0]
    info_P[ind_nii6548+2].limits=   [0., 0.]

    ; [NII] 6583
    ; redshift and width tied to narrow Ha 6563, flux tied to [NII] 6548
    info_P[ind_nii6583].tied    =   'P['+ $
                                    strcompress(string(ind_ha6563_n), /remove_all)+ $
                                    '] - alog(6564.61/6585.28)'
    info_P[ind_nii6583+1].tied  =   'P['+ $
                                    strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                    ']'
    info_P[ind_nii6583+2].tied =   '2.96 * P['+ $
                                    strcompress(string(ind_nii6548+2), /remove_all)+ $
                                    ']'

	    if not (keyword_set(freesii) and keyword_set(fixed2sii)) then begin
        info_P[ind_sii6716].tied    =   'P['+ $
                                        strcompress(string(ind_ha6563_n), /remove_all)+ $
                                        '] - alog(6564.61/6718.29)'
        info_P[ind_sii6716+1].tied  =   'P['+ $
                                        strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                        ']'
        info_P[ind_sii6716+2].limited=  [1, 0]
        info_P[ind_sii6716+2].limits=   [0., 0.]
;
;    ;; [SII] 6731
;    ;; redshift tied to narrow Ha 6563 and width tied to narrow-line model or be free
        info_P[ind_sii6731].tied    =   'P['+ $
                                        strcompress(string(ind_ha6563_n), /remove_all)+ $
                                        '] - alog(6564.61/6732.67)'
        info_P[ind_sii6731+1].tied  =   'P['+ $
                                        strcompress(string(ind_nlinemodel+1), /remove_all)+ $
                                        ']'

        info_P[ind_sii6731+2].limited   =  [1, 0]
        info_P[ind_sii6731+2].limits    =   [0., 0.]
    endif else begin

        info_P[ind_sii6731].limited    =   [1, 1]
        info_P[ind_sii6731].limits     =   [alog(6725.), alog(6740.)]
        info_P[ind_sii6731+1].limited  =   [1, 1]
        info_P[ind_sii6731+1].limits   =   [0.,500.]/(3d5*2.354)
        info_P[ind_sii6731+2].limited  =   [1, 0]
        info_P[ind_sii6731+2].limits   =   [0., 0.]

        info_P[ind_sii6716].tied    =   'P['+ $
                                        strcompress(string(ind_sii6731), /remove_all)+ $
                                        '] - alog(6732.67/6718.29)'
        info_P[ind_sii6716+1].tied  =   'P['+ $
                                        strcompress(string(ind_sii6731+1), /remove_all)+ $
                                        ']'
    endelse

	   ;;;;; broad components begin ;;;;;

    ; broad Ha 6563
    ; line center, width and flux constrainted
    info_P[ind_ha6563_b1].limited    =   [1, 1]
    ;info_P[ind_ha6563_b1].limits     =   [alog(6453), alog(6593)]
    info_P[ind_ha6563_b1].limits     =   [alog(6383), alog(6693)]
    info_P[ind_ha6563_b1+1].limited  =   [1, 1]
    info_P[ind_ha6563_b1+1].limits   =   [1000, 1d4]/(3d5*2.354)
    info_P[ind_ha6563_b1+2].limited  =   [1, 0]
    info_P[ind_ha6563_b1+2].limits   =   [0., 0.]

    info_P[ind_ha6563_b2].limited    =   [1, 1]
    ;info_P[ind_ha6563_b2].limits     =   [alog(6453), alog(6593)]
    info_P[ind_ha6563_b2].limits     =   [alog(6383), alog(6693)]
    info_P[ind_ha6563_b2+1].limited  =   [1, 1]
    info_P[ind_ha6563_b2+1].limits   =   [400, 1d4]/(3d5*2.354)
    info_P[ind_ha6563_b2+2].limited  =   [1, 0]
    info_P[ind_ha6563_b2+2].limits   =   [0., 0.]

	info_P[ind_hb4861_b1].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b1), /remove_all)+ $
                                    '] - alog(6564.61/4862.68)'
    info_P[ind_hb4861_b1+1].tied =   'P['+ $
                                    strcompress(string(ind_ha6563_b1+1), /remove_all)+ $
                                    ']'

    info_P[ind_hb4861_b1+2].limited   = [1, 0]
    info_P[ind_hb4861_b1+2].limits    = [0., 0.]

    info_P[ind_hb4861_b2].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b2), /remove_all)+ $
                                    '] - alog(6564.61/4862.68)'
    info_P[ind_hb4861_b2+1].tied =   'P['+ $
                                    strcompress(string(ind_ha6563_b2+1), /remove_all)+ $
                                    ']'
    info_P[ind_hb4861_b2+2].limited= [1, 0]
    info_P[ind_hb4861_b2+2].limits=  [0., 0.]
    ; broad H_gamma 4340
    ; redshift and width tied to broad Hbeta 4861
    info_P[ind_hc4340_b1].tied       =   'P['+ $
                                    strcompress(string(ind_ha6563_b1), /remove_all)+ $
                                    '] - alog(6564.61/4341.68)'
    info_P[ind_hc4340_b1+1].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b1+1), /remove_all)+ $
                                    ']'
    info_P[ind_hc4340_b1+2].limited  =   [1, 0]
    info_P[ind_hc4340_b1+2].limits   =   [0., 0.]

    info_P[ind_hc4340_b2].tied       =   'P['+ $
                                    strcompress(string(ind_ha6563_b2), /remove_all)+ $
                                    '] - alog(6564.61/4341.68)'
    info_P[ind_hc4340_b2+1].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b2+1), /remove_all)+ $
                                    ']'
    info_P[ind_hc4340_b2+2].limited  =   [1, 0]
    info_P[ind_hc4340_b2+2].limits   =   [0., 0.]

	; broad HeII 4686
    ; redshift and width tied to broad Hbeta 4861
    info_P[ind_heii4686_b].tied     =   'P['+ $
                                    strcompress(string(ind_ha6563_b1), /remove_all)+ $
                                    '] - alog(6564.61/4687.02)'
    info_P[ind_heii4686_b+1].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b1+1), /remove_all)+ $
                                    ']'
    info_P[ind_heii4686_b+2].limited  =   [1, 0]
    info_P[ind_heii4686_b+2].limits   =   [0., 0.]

    ; broad HeI 5876
    ; redshift and width tied to broad Halpha 6563
    info_P[ind_hei5876_b].tied     =   'P['+ $
                                    strcompress(string(ind_ha6563_b1), /remove_all)+ $
                                    '] - alog(6564.61/5877.29)'
    info_P[ind_hei5876_b+1].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b1+1), /remove_all)+ $
                                    ']'
    info_P[ind_hei5876_b+2].limited  =   [1, 0]
    info_P[ind_hei5876_b+2].limits   =   [0., 0.]

    ; broad FeII lines
    ; line width tied to broad Ha 6563
    info_P[ind_bl_FeII+1].tied   =   'P['+ $
                                    strcompress(string(ind_ha6563_b1+1), /remove_all)+ $
                                    '] /2.3548' ; FWHM of gauss is 2.3548 sigma
    info_P[ind_bl_FeII+2].limited  =   [1, 0]
    info_P[ind_bl_FeII+2].limits   =   [0., 0.]

	; broken powerlaw continuum

	    ind_fit = where(good_wave gt 4200. and good_wave lt 7200.) ;hyliu modified it to 4200
    if (ind_fit[0] eq -1) then begin
        print, 'no good points!'
        goto, skip
        ;stop
    endif
    if (n_elements(ind_fit) lt 100) then begin
        print, 'Good points less than 100!'
        goto, skip
        ;stop
    endif

    fit_x   = good_lnwave[ind_fit]
    fit_y   = good_flux[ind_fit]
    fit_err = good_error[ind_fit]
    ;fit_x   = good_lnwave
    ;fit_y   = good_flux
    ;fit_err = good_error

    time0   = systime(1)
    params  = mpfitexpr(model_fun, fit_x, fit_y, fit_err, $
                        Parinfo=info_P, $
                        yfit=yfit, perror=perror, $
                        bestnorm=bestnorm, dof=dof, $
                        status=status,maxiter=500,gtol=1.d-8,$
                        ftol=1.d-8,/quiet)
    time1   = systime(1)

	print, ''
    print, 'Fitting Time Consumed: ', time1-time0, 'Seconds', $
            format='(A25, D20.3, A8)'
    print, 'status:', format='(A, $)'
    print, status
    print, 'Chi:', format='(A, $)'
    ;print, bestnorm, dof, bestnorm/dof
    print, bestnorm
    print,"N: ",n_elements(fit_y),n_elements(good_flux[ind_fit])
