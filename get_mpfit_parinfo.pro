;+
; NAME:
;   get_mpfit_parinfo
;
; AUTHOR:
;   He-Yang Liu
;   My Github homepage:
;       https://github.com/AstroHeyang
;
; PURPOSE:
;	Get the mpfit expression for qso_fit. 
;
; CALLING SEQUENCE:
;   par_info = get_mpfit_parinfo()
;
; DESCRIPTION:
;   NULL
;
; INPUT:
;	line_infos    
;
; KEYWORD PARAMETERS:
;   NULL
;
; OUTPUT:
;   return a struct array 
;-

function get_mpfit_parinfo, line_infos 

	if not keyword_set(line_infos) then begin
		line_infos = get_line_info()
	endif

	nComponents = n_elements(line_infos)
	num_par = line_infos[nComponents-1].mpfit_expr_index + $
		line_infos[nComponents-1].npar
	par_info_str = {			$
		value:		0.,			$
        limited:    [0, 0],		$
        limits:     [0., 0.],	$
        tied:       ''			$
		}
    par_info = replicate(par_info_str, num_par)
	
	;; setting for power-law
	;; default: limited = [0, 0]
	npar0 = line_infos[0].npar
	par_info[0:npar-1].values = line_infos[0].values
	par_info[0].limited = [1, 0]
	par_info[0].limits = [0, 0]
	par_info[1:npar-1].limited = [1, 1]
	par_info[1:npar-1].limits = [-5, 5]

	;; setting for lines
	for i = 1l, nComponents-1 do begin
		npar = line_infos[i].npar
		this_index = line_infos[i].mpfit_expr_index

		;; set initial values
		;; TODO: change the initial values for nGauss model
		this_values = line_infos[i].values
        this_tied_single = line_infos[i].tied
		count_par = 0
		while count_par lt npar do begin
            par_info[this_index:this_index+2].value = this_values		
			this_tied = [this_tied, this_tied_single]
			this_index += 3
			count_par += 3
		endwhile

		;; set tied, limits and limited
		;; if this line is not tied, then limited and limits are set
        for j = 0, n_elements(this_tied) do begin
			if this_tied[j] eq -1 then begin
                line_index = line_infos[this_tied[j]]
                tied_index = line_infos[line_index].mpfit_expr_index
				par_info[this_index].tied = $
					'P[' + strcompress(string(tied_index+j), /remove_all) + $
					']'
			endif else begin
				

		endfor
		

	endfor
			
	for i = 0l, nComponents-1 do begin



