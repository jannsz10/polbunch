		*! polbunch version date 20241209
		* Author: Martin Eckhoff Andresen
		* This program is part of the polbunch package.
		
		/*
		Add support for constant approx?
		bootstrap -1 => no standard errors reported
		
		
		REPORTING
		bunching eq: B, shift, marginal response, (elasticity)
		estimator 0:  	report estimates from unrestricted model - notransform just sums over bunching coefs
						+ bunching eq w/o std. errors (input: B,h0,bw)
		estimator 1: 	Report h0,h1,B (w/ transform), h0,bunch coefs (w/o transform) + bunching eq w/o std. errors
		estimator 2: 	Report h0,h1,B,delta (w/transform), h0 (-cons), bunch coefs, delta (w/o transform) + bunching eq w/o std errors
		estimator 3: 	Report h0,h1,B + bunching eq w/ std errors (w/ transform), h0 (-cons), bunch coefs (w/o transform)
		*/
		cap prog drop polbunch
		program polbunch, eclass sortpreserve
			syntax varlist(min=1 max=2) [if] [in],  CUToff(real) [bw(numlist min=1 max=1 >0)  ///
			LIMits(numlist max=2 min=2 integer) ///
			t0(numlist min=1 max=1  >=0 <=1) ///
			t1(numlist min=1 max=1  >=0 <=1) ///
			POLynomial(integer 7) ///
			NOIsily ///
			estimator(integer 3) /// Specify estimator - 3 = theoretically consistent efficient estimator, 2 = chetty, 1 = no adjustment, 0=data to the left only
			nodrop ///
			bootstrap(integer 200) ///
			initvals(string) ///
			notransform ///
			nopositiveshift ///
			bootreps(integer 0) ///
			log ///
			]
			
			quietly {
				if !inlist(`estimator',0,1,2,3) {
						noi di as error "Option estimator can take only values 1 (using data to the left only), 2 (no adjustment), 2 (Chetty et. al. adjustment) or 3 (theoretically consistent and efficient estimator)."
						exit 301
				}
				tempvar touse
				marksample touse
				preserve
				drop if !`touse'

				if `polynomial'<=0 {
					noi di as error "Polynomial must be a positive integer"
					exit 301
				}
				
				//limits option
				if "`limits'"=="" {
					loc L=1
					loc H=0
				}
				else gettoken L H: limits
				
				loc zH=`cutoff'+`H'*`bw'
				loc zL=`cutoff'-`L'*`bw'
				
				// check varlist vs bw opts
				loc nvars: word count `varlist' 
				if (`nvars'==1&"`bw'"=="")|(`nvars'==2&"`bw'"!="") {
					noi di as error "Varlist must either contain 1 variable (earnings z) and option bw be specified (individual level data) or 2 variables (frequency y and earnings bin z) and option bandwith not be specified (pre-binned data)."
					exit 301
				}
				if `nvars'==2 { //find bw in pre-binned data
					loc y: word 1 of `varlist'
					loc z: word 2 of `varlist'
					
					sort `z'
					tempname tmp
					gen `tmp'=`z'-`z'[_n-1]
					su `tmp'
					if r(Var)>0.01*r(mean) {
						noi di as error "Bandwidth differs in pre-binned data"
						exit 301
					}
					else loc bw=r(mean)
					sum `y'
					loc N=r(sum)
				}
				else { //collapse data
					loc z `varlist'
					tempvar bin y
					
					if "`drop'"!="nodrop" {
						su `z'
						if abs(`cutoff'-floor((`z'-r(min))/`bw')*`bw'-r(min))>`bw'/10 {
							drop if `z'<`cutoff'-floor((`cutoff'-r(min))/`bw')*`bw'
						}
						if abs(`cutoff'+floor((r(max)-`cutoff')/`bw')*`bw'-r(max))>`bw'/10 {
							drop if `z'>`cutoff'+floor((r(max)-`cutoff')/`bw')*`bw'
						}
					}
					count
					loc N=r(N)
					gen `bin'=ceil((`z'-`cutoff')/`bw')*`bw'+`cutoff'-`bw'/2
					
					collapse (count) `y'=`z', by(`bin')
					rename `bin' `z'
				}
				
				//Put bins in e(table)
				tempname table
				mkmat `y' `z', matrix(`table')
				mat colnames `table'= freq `z'
				
				//check if there are people on either side
				count if `z'<`cutoff'
				if r(N)==0 {
					noi di as error "No individuals in sample allocates below cutoff."
					exit 301
					}
				
				count if `z'>`cutoff'
				if r(N)==0 {
					noi di as error "No individuals in sample allocates above cutoff."
					exit 301
					}
					
				//NAMES
				forvalues i=1/`polynomial' { 
					if "`rhsvars'"=="" loc rhsvars  c.`z'
					else loc rhsvars `rhsvars'##c.`z'
					loc coleq0 `coleq0' h0
					loc coleq1 `coleq1' h1
					}
			
				loc coleq0 `coleq0' h0
				loc coleq1 `coleq1' h1
				
				fvexpand `rhsvars'
				loc names `r(varlist)' _cons
				
				//gen dummies
				tempvar fw dupe dum dum2 bunch
				gen `dum'=`z'>`cutoff'
				gen `dum2'=`dum'
				egen `bunch'=group(`z') if inrange(`z',`zL',`zH')
				replace `bunch'=0 if `bunch'==.
				count
				loc numbins=r(N)
				
				//Evaluate OVB
				loc nmiss=1
				while `nmiss'>0 {
					_rmcoll `y' 0.`dum'#(`rhsvars') 1.`dum'#(`rhsvars') 1.`dum' b0.`bunch'
					loc nmiss=strpos("`=r(varlist)'","o.")>0
					if `nmiss'>0 {
						loc note note
						loc polynomial=`polynomial'-1
						if `polynomial'==0 {
							noi di in red "Omitted variable problems."
							exit 301
							}
						
						//NEW NAMES
						loc coleq0 
						loc coleq1 
						forvalues i=1/`polynomial' { 
							if `i'==1 loc rhsvars  c.`z'
							else loc rhsvars `rhsvars'##c.`z'
							loc coleq0 `coleq0' h0
							loc coleq1 `coleq1' h1
						}
				
						loc coleq0 `coleq0' h0
						loc coleq1 `coleq1' h1
						
						fvexpand `rhsvars'
						loc names `r(varlist)' _cons
						}
					}
				if "`note'"=="note" {
					noi di as text "Note: Polynomial order lowered to `polynomial' because of multicollinearity problems with the specified polynomial."
				}
					
				//ESTIMATE UNRESTRICTED MODEL (as benchmark or main model if estimator==0)
				`noisily' reg `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' b0.`bunch', nocons
				loc df_m_u=e(df_m)
				tempvar pred rss
				predict double `pred'
				gen double `rss'=`y'*(`N'-`pred')^2+(`N'-`y')*(-`pred')^2
				su `rss'
				loc rss_u=r(sum)
				
				//IF  ESTIMATOR==0, this is the results, just organize
				if `estimator'==0 {
					if `bootreps'==0 {
						varcorrect `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' b0.`bunch', smallsample
						mat `V'=r(V)
						if "`transform'"!="notransform" {
							ereturn repost V=`V'
							bunchcalc, estimator(0) polynomial(`polynomial') bw(`bw') H(`H') L(`L') b0(_b[b0]) t0(`t0') t1(`t1') `constant' `positiveshift' `log'
							mat `b'=r(b)
							mat `V'=r`V)
						}
						else mat `b'=e(b)
					}
					else { //bootstrap standard errors for estimator 0
						tempname p
						gen double `p'=`y'/`N'
						noi di as text "Performing bootstrap repetitions..."
						bunchcalc, estimator(0) polynomial(`polynomial') bw(`bw') H(`H') L(`L') b0(_b[b0]) t0(`t0') t1(`t1') `constant' `positiveshift' `log' boot
						mat `b'=r(b)
						loc nlcom r(nlcom)
						simulate, reps(`bootreps'): bssim, modstr(0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' b0.`bunch') p(`p') obs(`N') nlcom(`nlcom') `log' `positiveshift' `constant'
						
					}
				}
				
				else { //estimate with NLS if estimator>0
				//FIND GOOD STARTING VALUES FROM THE UNRESTRICTED MODEL		
				tokenize `names'
				
				if "`log'"=="" {
					loc shiftinit=(_b[1.`dum'#`z']/_b[0.`dum'#`z'])^(1/2)-1
					if "`positiveshift'"=="nopositiveshift" {
						if !mi(`shiftinit') loc init shift `shiftinit'
						else loc init shift 0.1
					}
					else {
						if !mi(`shiftinit')&`shiftinit'>0 loc init lnshift `=ln(`shiftinit')'
						else loc init lnshift `=ln(0.1)'
					}
				}
					else { //initvals if log - improve!!
						if "`positiveshift'"=="nopositiveshift" loc init shift 0.1
						else loc init lnshift `=ln(0.1)'
					}
				
				****** specify NL model string *****
				if "`positiveshift'"!="nopositiveshift" loc shift exp({lnshift})
				else loc shift {shift}
	
				// b0 parameter - integration constraint
				loc B=0
				forvalues b=1/`=`H'+`L'' {
					loc B `B'+{bunch`b'}
					}
				loc cns (`B')*`bw'
				if `estimator'!=2 {
					forvalues k=1/`polynomial' {
						if "`log'"=="" loc cns `cns' - (((`cutoff'*(1+`shift'))^(`=`k'+1')-`cutoff'^`=`k'+1')*({b`k'}/`=`k'+1'))
						else loc cns `cns' - (((`cutoff'+ln(1+`shift'))^(`=`k'+1')-`cutoff'^`=`k'+1')*({b`k'}/`=`k'+1'))
						}
					if "`log'"=="" loc b0 (`cns')/(`cutoff'*`shift')
					else loc b0 (`cns')/ln(1+`shift')
				}
				else {
					loc cns (`cns'/`shift')
					su `z'
					loc zmax=r(max)+`bw'/2				
					forvalues k=1/`polynomial' {
						loc cns `cns' - ({b`k'}/`=`k'+1')*(`zmax'^`=`k'+1'-`cutoff'^`=`k'+1')
					}
				loc b0 (`cns')/(`zmax'-`cutoff')
				}
				
				if `estimator'==0 loc modstr `b0'*(`z'<=`cutoff')
				else if inlist(`estimator',1,2)|(`estimator'==3&"`log'"=="log") loc modstr `b0'
				else loc modstr (`b0')*(1+`shift')^(`z'>`cutoff')
		
				/// b1-bK parameters, gk parameters for estimator 0
				tokenize `names'
				if `estimator'==0 {
					forvalues k=1/`polynomial' {
						loc modstr `modstr' + {b`k'}*(`z'<=`cutoff')*`z'^`k'+{g`k'}*(`z'>`cutoff')*`z'^`k'
						loc init `init' b`k' `=_b[0.`dum'#``k'']' g`k' `=_b[1.`dum'#``k'']'
					}
					loc modstr `modstr' +{g0}*(`z'>`cutoff')
					loc init `init' g0 `=_b[1.`dum']'
				}
				else if inlist(`estimator',1,2) {
					forvalues k=1/`polynomial' {
						loc modstr `modstr' + {b`k'}*`z'^`k'
						loc init `init' b`k' `=_b[0.`dum'#``k'']'
						}
				}
				else if `estimator' == 3 {
					forvalues k=1/`polynomial' {
						if "`log'"=="" loc modstr `modstr' +{b`k'}*(`z'^`k') * (1+`shift')^((`z'>`cutoff')*`=`k'+1')
						else loc modstr `modstr' + {b`k'}*((`z'+(`z'>`cutoff')*ln(1+`shift'))^`k')
						loc init `init' b`k' `=_b[0.`dum'#``k'']'
						}
				}
				

				if `estimator'==2 loc modstr (1/(1+exp({lnshift})*(`z'>`cutoff')))*(`modstr')
				
				// bunch  region dummies
				forvalues b=1/`=`H'+`L'' {
					loc modstr `modstr' + {bunch`b'}*`b'.`bunch'
					loc init `init' bunch`b' `=_b[`b'.`bunch']'
					}
					

				`noisily' nl (`y' = `modstr'), init(`init')
				
				tempname b V
				if "`transform'"!="notransform" {
					bunchcalc, estimator(`estimator') polynomial(`polynomial') bw(`bw') H(`H') L(`L') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log'
					mat `b'=r(b)
					if `bootreps'==0 mat `V'=r(V)
					else {
						loc nlcom r(nlcom)
						bssim ...
					}
				}
				estat ic
				loc aic=r(S)[1,5]
				loc df_m_r=e(df_m)
				tempname b V
				mat `b'=e(b)
				

				
				if `estimator'>0 {
					loc F=((`rss_r'-`rss_u')/(`df_m_u'-`df_m_r'))/(`rss_u'/(`N'-1))		
					loc p=Ftail(`df_m_u'-`df_m_r',`N'-1,`F')
				}
				
				if "`transform'"!="notransform" {
					loc names `names' `names' shift number_bunchers excess_mass marginal_response
					if "`t0'"!=""&"`t1'"!="" loc names `names' elasticity
				}
				else {
					tokenize `names'
					loc names
					forvalues i=1/`=`L'' {
						loc names `names' number_bunchers:`=`L'-`i'+1'.below
					}
					if `H'>0 {
						forvalues i=1/`=`H'' {
							loc names `names' number_bunchers:`i'.above
						}
					}
					if "`positiveshift'"!="nopositiveshift" loc names `names' shift:lnshift
					else loc names `names' shift:shift
					forvalues k=1/`polynomial' {
						loc names `names' h0:``k''
					}
					if `estimator'==0 {
						forvalues k=1/`polynomial' {
							loc names `names' h1:``k''
						}
						loc names `names' h1:_cons
					}
				}
				}
				
				mat colnames `b'=`names'
				mat colnames `V'=`names'
				mat rownames `V'=`names'
				if "`transform'"!="notransform" {
					mat coleq `b'=`coleq0' `coleq1' bunching
					mat coleq `V'=`coleq0' `coleq1' bunching
					mat roweq `V'=`coleq0' `coleq1' bunching
				}
				
				restore
				
				//return results
				eret post `b' `V', esample(`touse') obs(`N')
				if `estimator'>0 {
					estadd scalar F_mod=`F'
					estadd scalar p_mod=`p'
					estadd scalar F_df1=`=`polynomial'+1'
					estadd scalar F_df2=`N'-1
				}
				ereturn scalar polynomial=`polynomial'
				ereturn scalar bandwidth=`bw'
				ereturn scalar cutoff=`cutoff'
				ereturn scalar lower_limit=`zL'
				ereturn scalar upper_limit=`zH'
				ereturn scalar aic=`aic'
				ereturn local cmdname "polbunch"
				ereturn local title 	"Polynomial bunching estimates"
				ereturn local cmdline 	"polbunch `0'"
				ereturn matrix table=`table'
				ereturn local binname "`z'"
				ereturn scalar bw=`bw'
				if `bootreps'>0 estadd local vcetype "bootstrap"
				else estadd local vcetype "cluster"
				if "`log'"=="log" ereturn scalar log=1
				else ereturn scalar log=0
				/*foreach s in `scalarlist' {
					ereturn local `s'=``s''
				}*/
					
				//Display results
				noi {
					di _newline
					di "`e(title)'"
					eret di
					if `estimator'>0 {
						di "F-test of model assumptions: {col 42}F(`=`df_m_u'-`df_m_r'',`=`N'-1') test statistic {col 72}`: di %12.4f `F''"
						di "{col 42}p-value {col 72}`: di %12.4f `p''"
						di "{hline 83}"
						}
					if inlist(`estimator',1,2) {
						di "Note: Estimator is not consistent with iso-elastic labor supply model and estimates are biased."
					}
					}
			}
				

		end
		
program define varcorrect, rclass
syntax anything, [nl smallsample]
	
	qui {
		gettoken y anything: anything
		tempvar res rss
		tempname g V
		loc numbins=_N
		su `y'
		loc N=r(sum)
		loc k=e(df_m)
		//if "`nl'"=="" loc ++k
		if "`smallsample'"=="" loc smallsample=1
		else loc smallsample=sqrt((`N'/(`N'-1))*((`N'*`numbins'-1)/(`N'*`numbins'-`k')))
		if "`nl'"=="nl" {
			predictnl `res'=predict(), g(`g') iterate(1000)
			gen `rss'=`y'*(`N'-`res')^2+(`N'-`y')*(`res')^2
			replace `res'=-sqrt(`smallsample')*`res' //small sample adj
			mata: st_matrix("`V'",varcorrect(st_data(.,"`g'*"),st_data(.,"`y'"),st_data(.,"`res'"),0,`smallsample'))
			}
		else {
			predict `res'
			gen `rss'=`y'*(`N'-`res')^2+(`N'-`y')*(`res')^2
			replace `res'=-sqrt(`smallsample')*`res' //small sample adj
			mata: st_matrix("`V'",varcorrect(st_data(.,"`anything'"),st_data(.,"`y'"),st_data(.,"`res'"),0,`smallsample'))
		}
		
		return matrix V=`V'
		su `rss'
		return local rss=r(sum)
	}
	end

mata:
	
function varcorrect(real matrix X, real matrix fw, real matrix e, real scalar addcons, real scalar smallsample)
	{
	B=length(fw)
	N=sum(fw)
	if (addcons==1) X=X,J(rows(X),1,1)
	meat=0
	for (i=1; i<=B; i++) {
		e[i]=e[i]+N*smallsample
		meat=meat:+fw[i]:*(X' * e * e' * X)
		e[i]=e[i]-N*smallsample
	}
	return(invsym(quadcross(X,X):*N)*meat*invsym(quadcross(X,X):*N))
	}
				
function eresp(real scalar B,real scalar tau,real matrix cf, real scalar bw)
	{
	integral=polyinteg(cf,1)
	integral[1]=-polyeval(integral,tau)-B*bw
	roots=polyroots(integral)
	realroots=Re(select(roots, Im(roots):==0))
	out=sort(select(realroots,realroots:>tau)',1)'
	if (cols(out)==0) {
		return(.)
	}
	else return(out[1]-tau)
	}
			
end

// bunchcalc: Transforms estimates to h0,h1,B,excess_mass,shift,marginal_response,(elasticity)
program bunchcalc, rclass
	syntax, estimator(integer 3) polynomial(integer) bw(real) H(integer) L(integer) b0(string) [t0(real) t1(real) constant nopositiveshift boot log nlcom(string)]
	
	//compute nlcom string
	if "`nlcom'"=="" {
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
					else loc nlcom `nlcom' (ln(1+exp((`bw'*`B')/(`m'))/`cutoff')/(ln(1-`t0')-ln(1-`t1'))) 
				}
		}
		else if "`constant'"!="" { //calculcate shift/MR/el also if estimator <3 using the constant approx
			if "`log'"=="" loc nlcom `nlcom' (((`bw'*`B')/(`m'))/`cutoff')
			else loc nlcom `nlcom' (exp((`bw'*`B')/(`m')-`cutoff'))
			loc nlcom `nlcom' ((`bw'*`B')/(`m'))
			if "`t0'"!=""&"`t1'"!="" {
				if "`log'"==""  (ln(((`bw'*`B')/(`m'))/`cutoff'+1)/(ln(1-`t0')-ln(1-`t1')))
				else loc nlcom `nlcom' (ln(1+exp((`bw'*`B')/(`m')-`cutoff'))/(ln(1-`t0')-ln(1-`t1')))
			}
		}
	}

	/// TRANSFORM ESTIMATES
	tempname b
	if "`boot'"=="" {
		tempname V
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
				mat `b'=`b',`=exp(eresp-`cutoff')',eresp	
				if "`t0'"!=""&"`t1'"!="" mat `b'=`b',`=ln(exp(eresp-`cutoff')+1)/(ln(1-`t0')-ln(1-`t1'))'	
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
				mat `b'=`b',`=exp(eresp-`cutoff')',eresp		
				if "`t0'"!=""&"`t1'"!="" mat `b'=`b',`=ln(exp(eresp-`cutoff')+1)/(ln(1-`t0')-ln(1-`t1'))'	
			}
		}
	}
	
	//NAMES!!
	
	return local nlcom="`nlcom'"
	return matrix b=`b'
	if "`boot'"=="" return matrix V=`V'
	
	end

// Helper program for the binned bootstrap
program bssim, eclass
	syntax , estimator(integer) polynomial(integer) bw(real) H(integer) L(integer) b0(string) p(varname) obs(real) estopts(string) [t0(real) t1(real) nl nlcom(string) log constant nopositiveshift boot notransform]
	preserve
	tempvar y
	
	gen `y'=.
	loc i=0
	while `obs'>0&`i'<`=_N-1' {
		loc ++i
		replace `y'=max(0,`obs'-rbinomial(`obs',1-`p')) in `i'
		replace `p'=`p'/(1-`p'[`i']) if _n>`i'
		loc obs=`obs'-`y'[`i']
	}
	if `obs'>0 replace `y'=`obs' in `=_N'
	else replace `y'=0 if `y'==.
	
	if "`nl'"=="" reg `y' `modstr', `estopts'
	else nl (`y'=`modstr'), `estopts'
	
	tempname b
	if "`transform'"=="notransform" {
		bunchcalc, estimator(`estimator') polynomial(`polynomial') bw(`bw') H(`H') L(`L') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log' boot nlcom(`nlcom')
		mat `b'=r(b)
	}
	else mat `b'=e(b)
	ereturn post `b', obs(`obs')
	end
