program polbunchgendata, rclass
	syntax , obs(real) zstar(real) el(real) t0(real) t1(real) [distribution(string) log]
	
	quietly {
	clear
	set obs `obs'
	
	if !inrange(`t0',0,1)|!inrange(`t1',0,1)|`t1'<`t0'|`el'<0 {
		noi di as error "tax rates t1 and t0 must be between 0 and 1, t1>t0, and the elasticity must be nonnegative. "
	}
	tempname r
	if "`distribution'"=="" loc distribution triangular(0,3,0)
	if strpos("`distribution'","triangular")>0 {
		cap {
		loc a=substr("`distribution'",12,1)
		loc b=substr("`distribution'",14,1)
		loc c=substr("`distribution'",16,1)
		gen double `r'=runiform()
		gen double z=`a'+sqrt(`r'*(`b'-`a')*(`c'-`a')) if `r'<(`c'-`a')/(`b'-`a')
		replace z=`b'- sqrt((1-`r')*(`b'-`a')*(`b'-`c')) if `r'>(`c'-`a')/(`b'-`a')
		}
	}
	else cap gen double z=`distribution'
	if _rc!=0 {
		noi di as error "Use a valid stata random number distribution function (with parameters) in distribution(). Additionally, triangular(a,b,c) is allowed."
	}
	su z
	if !inrange(`zstar',r(min),r(max)) {
		noi di as error "Cutoff=`zstar' is not within the support specified by the distribution in distribution()."
		exit 301
	}
	if "`log'"=="log" {
		replace z=`zstar' if inrange(z,`zstar',`zstar'+`el'*(ln(1-`t0')-ln(1-`t1')))
		replace z=z-`el'*(ln(1-`t0')-ln(1-`t1')) if z>`zstar'
	}
	else {
		replace z=`zstar' if inrange(z,`zstar',`zstar'*((1-`t0')/(1-`t1'))^`el')
		replace z=z*((1-`t1')/(1-`t0'))^`el' if z>`zstar'
	}
	}
	
end
	