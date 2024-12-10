program polbunchsim, rclass
	syntax , obs(real)  zstar(real) bw(real) pol(real) el(real) t0(real) t1(real) [notransform distribution(string) log]
	qui {
	polbunchgendata , distribution(`distribution') obs(`obs') zstar(`zstar') el(`el') t0(`t0') t1(`t1') `log'	
	polbunch z, cutoff(`zstar') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `opts' `transform'
	if "`transform'"=="notransform" nlcom (elasticity: ln(1+exp(_b[shift:lnshift]))/(ln(1-`t0')-ln(1-`t1'))), post
	test elasticity=`el'
	return scalar p=r(p)
	}
end