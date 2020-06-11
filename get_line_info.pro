;+
; NAME:
;   get_line_info
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;       https://github.com/AstroHeyang
;
; PURPOSE:
;	Defining the basic information of the components for the model of qso_fit. 
;
; CALLING SEQUENCE:
;   line_infos = get_line_info()
;
; DESCRIPTION:
;   NULL
;
; INPUT:
;	n_lines   
;
; KEYWORD PARAMETERS:
;   NULL
;
; OUTPUT:
;   return a struct array including the default parameters for mpfit
;-
function get_line_info, par=par

	;; values and tied : [centroid, sigma, area] for Gaussian
	;; most of the limiteds are [1, 0], and the corresponding limits are [0, 0]
	;; the original struct is 
	;; line_info_struct = {
	;;	name: ' ', values: [0.d, 0.d, 0.d], tied: [-1, -1, -1], flux_ratio:1.0, is_broadline: 0, $
	;;	limited:[[0, 0], [0, 0], [0, 0]], limits:[[0, 0], [0, 0], [0, 0]] $
	;;	}
	;; this can be simplified to the following one.

	if not keyword_set(par) then begin
		par = read_par('default.par')
	endif

	n_lines = 25
	line_info_struct = { base, $
		name: '', values: [0.d, 0.d, 0.d], tied: [-1, -1, -1], flux_ratio: 1.0, is_broadline: 0, $
		limits_centroid:[0., 0.], model: 'g', mpfit_expr_index: 0, npar: 3 $
		}

	line_infos = replicate(line_info_struct, n_lines)
	line_infos[0:n_lines-1] = [ $
		{base, 'broken_powerlaw',	[1,		  -1.7, -1.7], [-1, -1, -1],  1.0, 0, [1.0,     1.0],'g', 0, 3}, $ ;0
		{base, 'nl_FeII_func'	,	[5019.83,  500,   10], [-1, -1, -1],  1.0, 0, [5018., 5020.],'g', 0, 3}, $ ;1
		{base, 'bl_FeII_func'	,	[5019.83, 3000,   10], [-1, -1, -1],  1.0, 1, [5018., 5020.],'g', 0, 3}, $ ;2 
		{base, 'nl_Hgamma_4340'	,	[4341.68,  300,  200], [11, 11, -1],  1.0, 0, [4259., 4449.],'g', 0, 3}, $ ;3
        {base, 'bl_Hgamma_4340'	,	[4341.68, 2000,  500], [20, 20, -1],  1.0, 1, [4259., 4449.],'g', 0, 3}, $ ;4
        {base, 'OIII_4363'		,	[4364.44,  300,   50], [11, 11, -1],  1.0, 0, [4259., 4449.],'g', 0, 3}, $ ;5
        {base, 'nl_HeII_4686'	,	[4687.02,  300,   50], [11, 11, -1],  1.0, 0, [4642., 4764.],'g', 0, 3}, $ ;6
        {base, 'bl_HeII_4686'	,	[4687.02, 2000,  500], [20, 20, -1],  1.0, 1, [4642., 4764.],'g', 0, 3}, $ ;7
        {base, 'nl_Hbeta_4861'	,	[4862.68,  300,  200], [11, 11, -1],  1.0, 0, [4744., 5135.],'g', 0, 3}, $ ;8
        {base, 'bl_Hbeta_4861'	,	[4862.68, 2000, 2000], [20, 20, -1],  1.0, 1, [4744., 5135.],'g', 0, 3}, $ ;9
        {base, 'OIII_4959'		,	[4960.30,  500,  100], [11, 11, 11], 0.33, 0, [4744., 5135.],'g', 0, 3}, $ ;10
        {base, 'OIII_5007'		,	[5008.24,  500,  300], [-1, -1, -1],  1.0, 0, [4744., 5135.],'g', 0, 3}, $ ;11
        {base, 'NI_5199'		,	[5199.35,  500,  300], [19, 19, -1],  1.0, 0, [5190., 5222.],'g', 0, 3}, $ ;12
        {base, 'nl_HeI_5876'	,	[5877.29,  300,   50], [19, 19, -1],  1.0, 0, [5741., 6010.],'g', 0, 3}, $ ;13
        {base, 'bl_HeI_5876'	,	[5877.29, 2000,  500], [20, 20, -1],  1.0, 1, [5741., 6010.],'g', 0, 3}, $ ;14
        {base, 'OI_6300'		,	[6302.05,  300,   50], [-1, -1, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}, $ ;15
        {base, 'SIII_6313'		,	[6312.06,  300,   50], [15, 15, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}, $ ;16
        {base, 'OI_6364'		,	[6365.54,  300,   50], [15, 15, 15], 0.33, 0, [6259., 6767.],'g', 0, 3}, $ ;17
        {base, 'NII_6548'		,	[6549.85,  300,   50], [19, 19, 21], 0.35, 0, [6259., 6767.],'g', 0, 3}, $ ;18
        {base, 'nl_Halpha_6563'	,	[6564.61,  500,  200], [-1, -1, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}, $ ;19
        {base, 'bl_Halpha_6563'	,	[6564.61, 2000, 2000], [-1, -1, -1],  1.0, 1, [6259., 6767.],'g', 0, 3}, $ ;20
        {base, 'NII_6583'		,	[6585.28,  300,  150], [19, 19, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}, $ ;21
        {base, 'HeI_6678'		,	[6680.0,   300,  150], [19, 19, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}, $ ;22
        {base, 'SII_6716'		,	[6718.29,  300,   50], [19, 19, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}, $ ;23
        {base, 'SII_6731'		,	[6732.67,  300,   50], [19, 19, -1],  1.0, 0, [6259., 6767.],'g', 0, 3}  $ ;24
		]

	line_infos[1:n_lines-1].values[0] = alog(line_infos[1:n_lines-1].values[0])
	line_infos[1:n_lines-1].values[1] = line_infos[1:n_lines-1].values[1] / (3D+5*2.354)
	line_infos[1:n_lines-1].values[2] = line_infos[1:n_lines-1].values[2] / line_infos[1:n_lines-1].values[0]
	line_infos.limits_centroid = alog(line_infos.limits_centroid)

	line_names = strupcase(line_infos.name)
    par_tags = tag_names(par)
    for i = 1, n_lines-1 do begin
		line_infos[i].mpfit_expr_index = line_infos[i-1].mpfit_expr_index + line_infos[i-1].npar
        index_this_tag = where(par_tags eq line_names[i])

        ;; current_index = last_index + nGauss * 3, nGauss is obtained from
        ;; the par file, else nGauss = 1 (default)
        if index_this_tag[0] ne -1 then begin
            line_infos[i].npar = par.(index_this_tag) * 3
        endif
    endfor

	;; Lorentz or not for broad lines
	broad_line_index = where(line_infos.is_broadline eq 1)
	if par.lorentz eq 1 then begin
		if broad_line_index[0] ne -1 then begin
			line_infos[broad_line_index].model = 'l'
		endif
	endif

	return, line_infos
end
