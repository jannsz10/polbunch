program polbunchsim, eclass
	 syntax [, log reps(integer 1) obs(integer 5000) cutoff(real 1) el(real 0.4) t0(real 0.2) t1(real 0.6) bw(real 0.01) bootreps(integer 500) pol(integer 1) notransform distribution(string) nopositiveshift estimator(numlist integer) btype(numlist integer) clist(string)]
	
	quietly {
	if "`btype'"=="" loc btype=1
	if "`estimator'"=="" loc estimator=3

	if "`clist'"=="" loc clist `"noconstant"'
	loc numc: word count `clist'
	if `numc'>2 {
	 	noi di in red "clist() can contain at most two strings"
		exit 301
	}
	tokenize `clist'
	forvalues i=1/`numc' {
	 	if !inlist("``i''","constant","noconstant") {
			noi di in red "clist() can contain only constant or noconstant"
			exit 301
		}
	}
	loc numb: word count `btype'
	loc nume: word count `estimator'
	loc numc: word count `clist'
	loc numest=`numb'*`nume'*`numc'
	tempname b V z
	polbunchgendata `z', obs(`obs') cutoff(`cutoff') el(`el') t0(`t0') t1(`t1') `log'
		
	foreach bt of numlist `btype'  {
		foreach e of numlist `estimator' {
			foreach c in `clist' {
			if "`c'"=="constant" loc cval=1
			else loc cval=0
			timer on 1
			if `bt'==0 polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' estimator(`e') bootreps(0) `c'
			else if `bt'==1 polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' estimator(`e') bootreps(1) `c'
			else if `bt'==2 bootstrap, reps(`bootreps'): polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' estimator(`e') bootreps(0) `c'
			else if `bt'==3 polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' `constant' `postiveshift' estimator(`e') bootreps(`bootreps') nobayes `c'
			else if `bt'==4 polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' estimator(`e') bootreps(`bootreps') nozero nobayes `c'
			else if `bt'==5 polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' estimator(`e') bootreps(`bootreps') `c'
			else if `bt'==6 polbunch `z', cutoff(`cutoff') pol(`pol') bw(`bw') t0(`t0') t1(`t1') `log' estimator(`e') bootreps(`bootreps') nozero `c'
			timer off 1
			timer list
			timer clear
			loc time=r(t1)
			cap confirm matrix e(V) 
			loc novar=_rc
			if `novar'==0 {
				cap confirm scalar _se[bunching:elasticity]
				loc novar=_rc
				if `novar'==0 {
					test elasticity=`el'
					loc p=r(p)
				}
			}
		
			//combine parameters of interest in b if running multiple estimation methods
			if `numest'>1 {
				mat `b'=nullmat(`b'),_b[bunching:elasticity],`time'
				loc names `names' se_`bt'_`e'_`cval' t_`bt'_`e'_`cval'
				if `novar'==0 {
					cap confirm scalar _se[bunching:elasticity]
					mat `b'=nullmat(`b'),_se[bunching:elasticity],`p'
					loc names `names' se_`bt'_`e'_`cval' p_`bt'_`e'
				}
				cap confirm scalar e(p_mod)
				if _rc==0 {
					mat `b'=nullmat(`b'),e(p_mod)
					loc names `names' p_mod_`bt'_`e'_`cval'
					}

				}
			}
		}
	}
	if `numest'>1 {
		mat colnames `b'=`names'
		ereturn post `b'
	}
	else {	
		/*mat `b'=e(b)
		if e(p_mod)!=. loc p_mod=e(p_mod)
		loc  polynomial=e(polynomial)
		tempname tab
		mat `tab'=e(table)
		
		cap confirm matrix e(V) 
		if _rc==0 {
			mat `V'=e(V)
			ereturn post `b' `V'
		}
		else ereturn post `b'
		if "`p'"!="" ereturn scalar p=`p'
		if "`p_mod'"!="" ereturn scalar p_mod=`p_mod'
		ereturn matrix 
		ereturn scalar polynomial=`polynomial'
		ereturn scalar time=`time'
		
		ereturn local cmd "polbunch"
		ereturn scalar polynomial=e(polynomial)
		*/
		estadd scalar time=`time'
		noi eret di
	}
	drop `z'
	}
	
end