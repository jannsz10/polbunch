program polbunchsim, rclass
	syntax , obs(real)  zstar(real) bw(real) pol(real) el(real) t0(real) t1(real) [notransform distribution(string) log noisily nopositiveshift estimator(integer 3) noisily notest bootreps(integer 0)]
	 {
	polbunchgendata , distribution(`distribution') obs(`obs') zstar(`zstar') el(`el') t0(`t0') t1(`t1') `log'	
	`noisily' polbunch z, cutoff(`zstar') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `transform' `noisily' `positiveshift' estimator(`estimator') bootreps(`bootreps')
		return scalar e=_b[bunching:elasticity]
		if `bootreps'>0 {
			test elasticity=`el'
			return scalar p=r(p)
		}
	}
end