		*! polbunch version date 20251701
		* Author: Martin Eckhoff Andresen
		* This program is part of the polbunch package.
		
		cap prog drop polbunch
		program polbunch, eclass sortpreserve
			syntax varlist(min=1 max=2) [if] [in],  CUToff(real) [bw(numlist min=1 max=1 >0)  ///
			LIMits(numlist max=2 min=2 integer) ///
			t0(numlist min=1 max=1  >=0 <=1) ///
			t1(numlist min=1 max=1  >=0 <=1) ///
			POLynomial(integer 7) ///
			NOIsily ///
			ESTimator(integer 3) /// Specify estimator - 3 = theoretically consistent efficient estimator, 2 = chetty, 1 = no adjustment, 0=data to the left only
			nodrop ///
			INITvals(string) ///
			notransform ///
			nopositiveshift ///
			BOOTreps(integer 1) ///
			log ///
			constant ///
			nodots /// suppress dots for bootstrap progress
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
				
				if `bootreps'<0 {
					noi di as error "Option bootreps can only take values 0 (no inference), 1 (analytic standard errors, the default) or a positive integer >0 (binned bootstrap)."
					exit 301
				}
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
						loc min=r(min)
						loc max=r(max)
						if abs(`cutoff'-floor((`z'-`min')/`bw')*`bw'-`min')>`bw'/10 {
							drop if `z'<`cutoff'-floor((`cutoff'-`min')/`bw')*`bw'
						}
						if abs(`cutoff'+floor((`max'-`cutoff')/`bw')*`bw'-`max')>`bw'/10 {
							drop if `z'>`cutoff'+floor((`max'-`cutoff')/`bw')*`bw'
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
				gen byte `dum'=`z'>`cutoff'
				gen byte `dum2'=`dum'
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
				forvalues bval=1/`=`H'+`L'' {
					loc bunchvars `bunchvars' `bval'.`bunch'
					}
				`noisily' reg `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' `bunchvars', nocons
				loc df_m_u=e(df_m)
				tempvar pred rss
				predict double `pred'
				gen double `rss'=`y'*(`N'-`pred')^2+(`N'-`y')*(-`pred')^2
				su `rss'
				loc rss_u=r(sum)
				
				if `bootreps'>1 {
					tempname p
					gen double `p'=`y'/`N'`'
				}

				if `estimator'!=0 { //estimate with NLS if estimator>0
				//FIND GOOD STARTING VALUES FROM THE UNRESTRICTED MODEL		
				tokenize `names'
				
				if "`log'"=="" {
					loc shiftinit=(_b[1.`dum2'#`z']/_b[0.`dum'#`z'])^(1/2)-1
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
				
				if inlist(`estimator',1,2)|(`estimator'==3&"`log'"=="log") loc modstr `b0'
				else loc modstr (`b0')*(1+`shift')^(`z'>`cutoff')
		
				/// b1-bK parameters, gk parameters for estimator 0
				tokenize `names'
				if `estimator'==0 {
					forvalues k=1/`polynomial' {
						loc modstr `modstr' + {h0:``k''}*(`z'<=`cutoff')*`z'^`k'+{g`k'}*(`z'>`cutoff')*`z'^`k'
						loc init `init' b`k' `=_b[0.`dum'#``k'']' g`k' `=_b[1.`dum2'#``k'']'
					}
					loc modstr `modstr' +{g0}*(`z'>`cutoff')
					loc init `init' g0 `=_b[1.`dum2']'
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
				
				if `estimator'==2 loc modstr (1/(1+`shift'*(`z'>`cutoff')))*(`modstr')
				
				// bunch  region dummies
				forvalues b=1/`=`H'+`L'' {
					loc modstr `modstr' + {bunch`b'}*`b'.`bunch'
					loc init `init' bunch`b' `=_b[`b'.`bunch']'
					}
					
				
				//ESTIMATE NL MODEL
				`noisily' nl (`y' = `modstr'), init(`init')
				estat ic
				loc aic=r(S)[1,5]
				loc df_m_r=e(df_m)
				
				tempname res rss
				predict double `res'
				gen `rss'=`y'*(`N'-`res')^2+(`N'-`y')*(`res')^2
				su `rss'
				loc rss_r=r(sum)
				
				loc F=((`rss_r'-`rss_u')/(`df_m_u'-`df_m_r'))/(`rss_u'/(`N'-1))		
				loc pval=Ftail(`df_m_u'-`df_m_r',`N'-1,`F')
				}
				else { //repost b V from OLS with new names to prepare for inference + transform
					estat ic
					loc aic=r(S)[1,5]
					loc b0 _b[/b0]
					if "`transform'"!="notransform" {
						foreach l in b g {
							foreach k of numlist 1/`polynomial' 0 {
								loc newnames `newnames' /`l'`k'
							}
						}
						forvalues bval=1/`=`H'+`L'' {
							loc newnames `newnames' /bunch`bval'
							}
						}
					}
				
				//INFERENCE: BINNED BOOTSTRAP OR DELTA METHOD BASED ON CORRECTED VARIANCE, transformed or untransformed
				tempname b V
				if "`transform'"!="notransform"&`estimator'==0&`bootreps'!=1 {
						mat `b'=e(b)
						mat colnames `b'=`newnames'	
						eret post `b'
					}
				if `bootreps'==0 { //no inference
					if "`transform'"!="notransform" {
						bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`H') l(`L') b0(`b0') t0(`t0') t1(`t1') boot `constant' `positiveshift' `log'
						mat `b'=r(b)
						if r(exit)==1 {
							noi di in red "Could not find solution to polynomial equation for the response of the marginal buncher in one or more bootstrap repetitions. You could consider trying the constant approximation using the option "constant", or the option "notransform" to report raw estimates and then manually convert those to objects of interest post-estimation."
								exit 301
							}
						}
					else mat `b'=e(b)
				}
				else if `bootreps'==1 { //analytic standard errors
					if `estimator'>0 varcorrect `y', smallsample nl
					else varcorrect `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' `bunchvars', smallsample
					mat `V'=r(V)
					if "`transform'"!="notransform" {
						if `estimator'==0 {
							mat `b'=e(b)
							mat colnames `b'=`newnames'	
							mat rownames `V'=`newnames'
							mat colnames `V'=`newnames'
							eret post `b' `V'
						}
						else ereturn repost V=`V'
						bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`H') l(`L') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log'
						if r(exit)==1 {
							noi di in red "Could not find solution to polynomial equation for the response of the marginal buncher in one or more bootstrap repetitions. You could consider trying the constant approximation using the option "constant", or the option "notransform" to report raw estimates and then manually convert those to objects of interest post-estimation."
								exit 301
							}
						mat `b'=r(b)
						mat `V'=r(V)
					}
					else mat `b'=e(b)
				}
				else { //binned bootstrap
					if "`transform'"!="notransform" {
						bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`H') l(`L') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log' boot
						mat `b'=r(b)
						loc nlcom `=r(nlcom)'
					}
					else mat `b'=e(b)
					tempname tmpb
					forvalues s=1/`bootreps' {
						if `s'==1&"`dots'"!="nodots" nois _dots 0, title("Performing bootstrap repetitions...") reps(`bootreps')
						loc i=0
						loc factor=0
						loc obs=`N'
						while `obs'>0&`i'<`=_N-1' {
							loc ++i
							replace `y'=rbinomial(`obs',`p'/(1-`factor')) in `i'
							loc factor=`factor'+`p'[`i']
							loc obs=`obs'-`y'[`i']
						}
						if `obs'>0 replace `y'=`obs' in `=_N'
						else replace `y'=0 if _n>`i'
						
						if `estimator'>0 nl (`y'=`modstr'), init(`init')
						else reg  `y' 0.`dum'#(`rhsvars') 0.`dum' 1.`dum2'#(`rhsvars') 1.`dum2' `bunchvars', nocons
						
						if "`transform'"!="notransform" {
							if `estimator'==0 {
								mat `tmpb'=e(b)
								mat colnames `tmpb'=`newnames'	
								eret post `tmpb'
							}
							bunchcalc, estimator(`estimator') polynomial(`polynomial') cutoff(`cutoff') bw(`bw') h(`H') l(`L') b0(`b0') t0(`t0') t1(`t1') `constant' `positiveshift' `log' boot nlcom(`nlcom')
							if r(exit)==1 {
								noi di in red "Could not find solution to polynomial equation for the response of the marginal buncher in one or more bootstrap repetitions. You could consider trying the constant approximation using the option "constant", or the option "notransform" to report raw estimates and then manually convert those to objects of interest post-estimation."
								exit 301
							}
							mat `V'=nullmat(`V') \ r(b)
						}
						else mat `V'=nullmat(`V') \ e(b)
						if "`dots'"!="nodots" noi _dots `s' 0
						}
						
					clear
					svmat `V'
					corr _all, cov
					mat `V'=r(C)
					}
				

				//NAMES
				if "`transform'"!="notransform" {
					if `estimator'!=2 loc names `names' `names' number_bunchers excess_mass shift marginal_response
					else loc names `names' `names' delta number_bunchers excess_mass shift marginal_response
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
					
					if `estimator'>0 {
						if `estimator'==2 {
							if "`positiveshift'"!="nopositiveshift" loc names `names' delta:lndelta
							else loc names `names' delta:delta
						}
						else {
							if "`positiveshift'"!="nopositiveshift" loc names `names' shift:lnshift
							else loc names `names' shift:shift
						}
					forvalues k=1/`polynomial' {
						loc names `names' h0:``k''
					}
					if `estimator'==0 {
						loc names `names' h0:_cons
						forvalues k=1/`polynomial' {
							loc names `names' h1:``k''
						}
						loc names `names' h1:_cons
					}
				}
				}
				
				mat colnames `b'=`names'
				if `bootreps'>=1 {
					mat colnames `V'=`names'
					mat rownames `V'=`names'
				}
				if "`transform'"!="notransform" {
					mat coleq `b'=`coleq0' `coleq1' bunching
					if `bootreps'>=1 {
						mat coleq `V'=`coleq0' `coleq1' bunching
						mat roweq `V'=`coleq0' `coleq1' bunching
					}
				}
				
				restore
				
				//return results
				noi mat li `b'
				if `bootreps'>=1 eret post `b' `V', esample(`touse') depname(freq) obs(`N')
				else eret post `b', esample(`touse') obs(`N') depname(freq)
				if `estimator'>0 {
					estadd scalar F_mod=`F'
					estadd scalar p_mod=`pval'
					estadd scalar F_df1=`=`df_m_u'-`df_m_r''
					estadd scalar F_df2=`N'-1
				}
				ereturn scalar polynomial=`polynomial'
				ereturn scalar bandwidth=`bw'
				ereturn scalar cutoff=`cutoff'
				ereturn scalar lower_limit=`zL'
				ereturn scalar upper_limit=`zH'
				ereturn scalar aic=`aic'
				ereturn local cmd "polbunch"
				ereturn local title 	"Polynomial bunching estimates"
				ereturn local cmdline 	"polbunch `0'"
				ereturn matrix table=`table'
				ereturn local binname "`z'"
				ereturn scalar bw=`bw'
				if `bootreps'>1 estadd local vcetype "bootstrap"
				if "`log'"=="log" ereturn scalar log=1
				else ereturn scalar log=0
	
				//Display results
				noi {
					di _newline
					di "`e(title)'"
					eret di
					if `estimator'>0 {
						di "F-test of model assumptions: {col 42}F(`=`df_m_u'-`df_m_r'',`=`N'-1') test statistic {col 72}`: di %12.4f `F''"
						di "{col 42}p-value {col 72}`: di %12.4f `pval''"
						di "{hline 83}"
						}
					if inlist(`estimator',1,2) {
						di "Note: Estimator is not consistent with iso-elastic labor supply model and thus biased."
					}
					if "`constant'"!="" {
						di "Note: Using the constant approximation to the density to calculate the elasticity may lead to bias."
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
			replace `res'=-sqrt(`smallsample')*`res' //small sample adj
			
			mata: st_matrix("`V'",varcorrect(st_data(.,"`g'*"),st_data(.,"`y'"),st_data(.,"`res'"),0,`smallsample'))

			}
		else {
			predict `res'
			replace `res'=-sqrt(`smallsample')*`res' //small sample adj
			mata: st_matrix("`V'",varcorrect(st_data(.,"`anything'"),st_data(.,"`y'"),st_data(.,"`res'"),0,`smallsample'))
		}
		
		return matrix V=`V'
		gen `rss'=`y'*(`N'-`res')^2+(`N'-`y')*(`res')^2
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

