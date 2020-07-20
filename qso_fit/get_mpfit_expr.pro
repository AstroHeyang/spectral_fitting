;+
; NAME:
;   get_mpfit_expr
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
;   mpfit_expr = get_mpfit_expr()
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
;   return a string 
;-

function get_mpfit_expr, line_infos 

	if not keyword_set(line_infos) then begin
		line_infos = get_line_info()
	endif

	nComponents = n_elements(line_infos)
    mpfit_expr = $
		'broken_powerlaw(X, [P[0:2]]) + ' +			$
    	'nlFeIIfunc(X, [P[3:5]], /extended) + ' +	$
        'blFeIIfunc(X, [P[6:8]], /extended) + ' 

    for i = 3, nComponents-1 do begin
		this_index = line_infos[i].mpfit_expr_index
		npar = line_infos[i].npar
        if i lt nComponents-1 then begin
            mpfit_expr += $
				'ngauss(X, [P[' + strtrim(string(this_index), 2) + ':' + $
				strtrim(string(this_index + npar - 1), 2) + ']]) + '
        endif else begin
            mpfit_expr += $
				'ngauss(X, [P[' + strtrim(string(this_index), 2) + ':' + $
				strtrim(string(this_index + npar - 1), 2) + ']])'
        endelse
    endfor

    return, mpfit_expr
end
