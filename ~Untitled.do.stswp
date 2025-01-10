transform

est 0: 
input: h0,h1,B1,B2,B3 ..., cutoff, bw
bootoutput: h0,h1,B,excess mass,MR,shift,el
nlcomoutputmatrix: MR,shift,el
nlcomoutputstring: h0,h1,delta,B,excess mass

est 1
input: h0,B1,B2,B3 ..., cutoff, bw
bootoutput: h0,h1,B,excess mass,MR,shift,el
nlcomoutputmatrix: MR,shift,el
nlcomoutputstring: h0,h1,delta,B,excess mass

est 2:
input: h0 net of constant,B1,B2,B3 ..., cutoff, bw
bootoutput: h0,h1,delta,B,excess mass,MR,shift,el
nlcomoutputmatrix: MR,shift,el
nlcomoutputstring: h0,h1,delta,B,excess mass

est 3:
input h0 net of constant,shift,B1,B2,B3 ..., cutoff, bw
bootoutput: h0,h1,B,excess mass,MR,shift,el
nlcomoutputmatrix: empty
nlcomoutputstring h0,h1,B,excess mass,MR,shift,el

program bunchcalc, rclass
	syntax varname, estimator(integer 3) polynomial(integer) bw(real) H(integer) L(integer) boot bunch(varname) dum(varname) z(varname) [t0(real) t1(real)] constant nopositiveshift b0(string)
	
	//make expression for b
	
	//bk (parameters in h0)
	forvalues k=1/`=`polynomial'' {
		loc nlcom `nlcom' (_b[/b`k'])
	}

	//b0 (constant in h0)
	if "`positiveshift'"!="nopositiveshift" loc b0=subinstr("`b0'","{lnshift}","_b[/lnshift]",.)
	else loc b0=subinstr("`b0'","{shift}","_b[/shift]",.)
	forvalues k=1/`=`polynomial'' {
		loc b0=subinstr("`b0'","{b`k'}","_b[/b`k']",.)
	}
	forvalues b=1/`=`H'+`L'' {
		loc b0=subinstr("`b0'","{bunch`b'}","_b[/bunch`b']",.)
	}
	loc nlcom `nlcom' (`b0')

	//gk (parameters in h1)
	forvalues k=1/`=`polynomial'' {
		if `estimator'==0 loc nlcom `nlcom' (_b[/g`k'])
		else if `estimator'==1 loc nlcom `nlcom' (_b[/b`k'])
		else if `estimator'==2 loc nlcom `nlcom' (_b[/b`k']/(1+`shift'))
		else if "`log'"=="" loc nlcom `nlcom' (_b[/b`k']*(1+`shift')^(`=`k'+1'))
		else {
				loc str _b[/b`k']
				if (`polynomial'>`k') {
					forvalues n=`=`k'+1'/`=`polynomial'' {
						loc str `str' +_b[/b`n']*comb(`n',`k')*ln(1+`shift')^(`n'-`k')
					}
				}
			loc nlcom `nlcom' (`str')
			}
		}
	
	//g0 (constant in h1)
	if `estimator'==0 loc nlcom `nlcom' (_b[/g0])
	else if `estimator'==1 loc nlcom `nlcom' (`b0')
	else if `estimator'==2 loc nlcom `nlcom' ((`b0')/(1+`shift'))
	else if "`log'"=="" loc nlcom `nlcom' ((`b0')*(1+`shift')) 
	else {
		loc str `b0'
		forvalues n=1/`=`polynomial'' {
				loc str `str' +_b[/b`n']*ln(1+`shift')^`n'
			}
			loc nlcom `nlcom' (`str')
	}

	//number bunchers
	loc B=0
	forvalues b=1/`=`H'+`L'' {
		loc B `B'+_b[/bunch`b']
		}
	loc nlcom `nlcom' (`B')
	loc m `b0'
	forvalues k=1/`polynomial' {
		loc m `m' + _b[/b`k']*`cutoff'^(`k')
	}
	loc nlcom `nlcom' ((`B')/(`m')) //excess mass
	
	if `estimator'==3  { //nlcom also  shift, MR, el	
	
		//shift
		if "`positiveshift'"!="nopositiveshift" loc nlcom `nlcom' (exp(_b[/lnshift]))
		else loc nlcom `nlcom' (_b[/shift])
		
		//response of marginal buncher
		if "`log'"!="log" loc nlcom `nlcom' (`shift'*`cutoff')
		else loc nlcom `nlcom' (ln((1+`shift')*`cutoff')-ln(`cutoff'))
						
		//elasticity
		if "`t0'"!=""&"`t1'"!="" {
			if "`constant'"=="" loc nlcom `nlcom' (ln(1+`shift')/(ln(1-`t0')-ln(1-`t1')))
			else { //calc elasticity using constant approx
				if "`log'"=="" loc nlcom `nlcom' (ln(((`bw'*`B')/(`m'))/`cutoff'+1)/(ln(1-`t0')-ln(1-`t1'))) 
				else loc nlcom `nlcom' ((((`bw'*`B')/(`m'))/`cutoff')/(ln(1-`t0')-ln(1-`t1'))) 
			}
	}
	else if "`constant'"!="" { //calculcate shift/MR/el also if estimator <3 using the constant approx
		if "`log'"=="" loc nlcom `nlcom' (((`bw'*`B')/(`m'))/`cutoff')
		else loc nlcom `nlcom' (exp(eresp+`cutoff')/exp(`cutoff')-1) 
		loc nlcom `nlcom' ((`bw'*`B')/(`m'))
		if "`t0'"!=""&"`t1'"!="" {
			if "`log'"==""  (ln(((`bw'*`B')/(`m'))/`cutoff'+1)/(ln(1-`t0')-ln(1-`t1')))
			else loc nlcom `nlcom' ((`bw'*`B')/(`m')/(ln(1-`t0')-ln(1-`t1')))
	}

	tempname b V
	if "`boot'"=="" {
		nlcom `nlcom'
		mat `b'=r(b)
		mat `V'=r(V)
		if `estimator'<3&"`constant'"=="" {
			//calculate shift/MR/el without variance and add to b,V
			tempname h0
			forvalues k=1/`polynomial' {
				mat `h0'=nullmat(`h0'),_b[/b`k']
			}
			mat `h0'=`h0',`b0'
			mata: st_numscalar("eresp",eresp(`=`B'',`cutoff',st_matrix("`h0'"),`bw'))
			if "`log'"=="" {
				mat `b'=`b',`=eresp/`cutoff'',eresp
				if "`t0'"!=""&"`t1'"!="" mat `b'=`b',`=ln(eresp/`cutoff'+1)/(ln(1-`t0')-ln(1-`t1'))'
			}
			else {
				mat `b'=`b',`=exp(eresp+`cutoff')/exp(`cutoff')-1',eresp		
				if "`t0'"!=""&"`t1'"!="" mat `b'=`b',`=eresp/(ln(1-`t0')-ln(1-`t1'))'	
			}
			if "`t0'"!=""&"`t1'"!="" loc extra=3
			else loc extra=2
			mat `V'=[`V', J(rowsof(`V'),`extra',0) \ J(`extra',colsof(`V'),0) , J(`extra',`extra',0)]
		}
	}
	else {
		 while "`nlcom'"!="" {
			gettoken use nlcom: nlcom, match(parns)
			mat `b'=nullmat(`b'),`=`use''
		}
		if `estimator'<3&"`constant'"=="" {
			tempname h0
			forvalues k=1/`polynomial' {
				mat `h0'=nullmat(`h0'),_b[/b`k']
			}
			mat `h0'=`h0',`b0'
			mata: st_numscalar("eresp",eresp(`=`B'',`cutoff',st_matrix("`h0'"),`bw'))
			if "`log'"=="" {
				mat `b'=`b',`=eresp/`cutoff'',eresp
				if "`t0'"!=""&"`t1'"!="" mat `b'=`b',`=ln(eresp/`cutoff'+1)/(ln(1-`t0')-ln(1-`t1'))'
			}
			else {
				mat `b'=`b',`=exp(eresp+`cutoff')/exp(`cutoff')-1',eresp		
				if "`t0'"!=""&"`t1'"!="" mat `b'=`b',`=eresp/(ln(1-`t0')-ln(1-`t1'))'	
			}
		}
	}
	
	//NAMES!!
	
	return matrix b=`b'
	if "`boot'"=="" return matrix V=`V'
	
	end