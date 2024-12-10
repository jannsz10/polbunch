

program polbunchsim, rclass
	syntax , obs(real) a(real) b(real) cval(real) zstar(real) bw(real) pol(real) el(real) t0(real) t1(real) notransform [log]
	qui {
	clear
	set obs `obs'
	tempname r
	gen double `r'=runiform()
	gen double z=`a'+sqrt(`r'*(`b'-`a')*(`cval'-`a')) if `r'<(`cval'-`a')/(`b'-`a')
	replace z=`b'- sqrt((1-`r')*(`b'-`a')*(`b'-`cval')) if `r'>(`cval'-`a')/(`b'-`a')
	if "`log'"=="log" {
		replace z=`zstar' if inrange(z,`zstar',`zstar'+`el'*(ln(1-`t0')-ln(1-`t1')))
		replace z=z-`el'*(ln(1-`t0')-ln(1-`t1')) if z>`zstar'
	}
	else {
		replace z=`zstar' if inrange(z,`zstar',`zstar'*((1-`t0')/(1-`t1'))^`el')
		replace z=z*((1-`t1')/(1-`t0'))^`el' if z>`zstar'
	}
	
	polbunch z, cutoff(`zstar') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `opts' `transform'
	if "`transform'"=="notransform" nlcom (elasticity: ln(1+exp(_b[shift:lnshift]))/(ln(1-`t0')-ln(1-`t1'))), post
	test elasticity=`el'
	return scalar p=r(p)
	}
end