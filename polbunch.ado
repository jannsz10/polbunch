		*! polbunch version date 20241209
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
			log ///
			nointconstraint ///
			fill ///
			nodrop ///
			bootstrap(integer 200) ///
			chetty ///
			initvals(string) ///
			notransform ///
			positiveshift ///
			]
			
			
			quietly {

				tempvar touse
				marksample touse
				preserve
				drop if !`touse'

				if `polynomial'<=0 {
					noi di as error "Polynomial must be a positive integer"
					exit 301
				}
				
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
					if r(Var)!=0 {
						noi di as error "Bandwidth differs in pre-binnd data"
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
					noi di "`cutoff'"
					noi two bar `y' `z', barwidth(`bw')
				}
				
				//Put bins in e(table)
				tempname table
				mkmat `y' `z', matrix(`table')
				mat colnames `table'= freq `z'
				
				//limits option
				if "`limits'"=="" {
					loc L=1
					loc H=0
				}
				else gettoken L H: limits
				
				loc zH=`cutoff'+`H'*`bw'
				loc zL=`cutoff'-`L'*`bw'
				
				
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
				tempvar fw dupe dum bunch
				gen `dum'=`z'>`cutoff'
				egen `bunch'=group(`z') if inrange(`z',`zL',`zH')
				replace `bunch'=0 if `bunch'==.
				
				//expand data & gen dummies
				count
				loc numbins=r(N)
				expand 2, gen(`dupe')
				gen `fw'=`y' if `dupe'==0
				replace `fw'=`N'-`y' if `dupe'==1
				replace `y'=`N' if `dupe'==0
				replace `y'=0 if `dupe'==1
					
				//ESTIMATE UNRESTRICTED MODEL AS benchmark (++lower polynomial of omitted vars)
				loc nmiss=1
				while `nmiss'>0 {
					`noisily' reg `y' 0.`dum'#(`rhsvars') 1.`dum'#(`rhsvars') 1.`dum' b0.`bunch' [fw=`fw']
					loc bnames: colnames e(b)
					loc nmiss=strpos("`bnames'","o.")
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
					
				varcorrect `y' 0.`dum'#(`rhsvars') 1.`dum'#(`rhsvars') 1.`dum' b0.`bunch', fw(`fw') numbins(`numbins')

				tempname b V
				mat `b'=e(b)
				mat `V'=r(V)
				
				loc bnames: colnames `b'
				tokenize `bnames'
				forvalues k=1/`=colsof(`b')-1' {
					test ``k'', accum
				}
				loc F=r(F)
				ereturn repost b=`b' V=`V'
				estadd scalar N_clust=`N', replace
				estadd local vce="cluster", replace
				estadd scalar df_r =  `N'-1, replace
				estadd scalar df_m = 0, replace	
				estadd matrix beta=., replace
				estadd scalar F=`F', replace
				est sto unres
				
				
				//FIND GOOD STARTING VALUES FROM UNRESTRICTED MODEL		
				tokenize `names'
				
				if "`log'"=="" {
					loc shiftinit=(_b[1.`dum'#`z']/_b[0.`dum'#`z'])^(1/2)-1
					if "`positiveshift'"=="" {
						if !mi(`shiftinit') loc init shift `shiftinit'
						else loc init shift 0.1
					}
					else {
						if !mi(`shiftinit')&`shiftinit'>0 loc init lnshift `=ln(`shiftinit')'
						else loc init lnshift `=ln(0.1)'
					}
				}
					else { //initvals if log - improve!!
						if "`positiveshift'"=="" loc init shift 0.1
						else loc init lnshift `=ln(0.1)'
					}
				
				//specify NL model string
				if "`positiveshift'"=="positiveshift" loc shift exp({lnshift})
				else loc shift {shift}
				if "`intconstraint'"=="nointconstraint" {
					loc modstr {b0}
					loc init `init' `=_b[_cons]'
					}
				else {
					loc zmaxadj `cutoff'*(1+`shift')
					loc B=0
					loc lower=`cutoff'
					forvalues b=1/`=`H'+`L'' {
						loc B `B'+{bunch`b'}
						}
					loc cns (`B')*`bw'
					forvalues k=1/`polynomial' {
						loc cns `cns' - ({b`k'}/`=`k'+1')*((`cutoff'*(1+`shift'))^(`=`k'+1')-(`cutoff'^`=`k'+1'))
					}
					loc b0 (`cns')/(`cutoff'*`shift')
					loc modstr `b0'
					//else loc modstr (`b0')*(1+`shift')^(`z'>`cutoff')
					
			}
		
				tokenize `names'
				forvalues k=1/`polynomial' {
					if "`log'"=="" loc modstr `modstr' +{b`k'}*(`z'^`k')*(1+`shift')^((`z'>`cutoff')*`=`k'+1')
					else loc modstr `modstr' +{b`k'}*(`z'+(`z'>`cutoff')*ln(1+`shift'))^`k'
					loc init `init' b`k' `=_b[0.`dum'#``k'']'
					}
				
				forvalues b=1/`=`H'+`L'' {
					loc modstr `modstr' + {bunch`b'}*`b'.`bunch'
					loc init `init' bunch`b' `=_b[`b'.`bunch']'
					}

			
				`noisily' nl (`y' = `modstr') [fw=`fw'], init(`init')
				
				varcorrect `y' test, fw(`fw') numbins(`numbins') nl
				tempname b V
				mat `b'=e(b)
				mat `V'=r(V)
				
				loc bnames: colnames `b'
				tokenize `bnames'
				forvalues k=1/`=colsof(`b')-1' {
					test ``k'', accum
				}
				loc F=r(F)
				ereturn repost b=`b' V=`V'
				estadd scalar N_clust=`N', replace
				estadd local vce="cluster", replace
				estadd scalar df_r =  `N'-1, replace
				estadd scalar df_m = 0, replace	
				estadd matrix beta=., replace
				estadd scalar F=`F', replace
				//noi eret di
				est sto res
				
				lrtest res unres, force				
				loc chi2=r(chi2)
				loc p=r(p)
				loc df_chi2=r(df)
				
				estat ic
				loc aic=r(S)[1,5]
				
				//Compute transformations of interest using nlcom

				if "`transform'"!="notransform" {
					if "`positiveshift'"!="" loc shift=subinstr("`shift'","{lnshift}","_b[/lnshift]",.)
					else loc shift=subinstr("`shift'","{shift}","_b[/shift]",.)
					forvalues k=1/`=`polynomial'' {
						loc nlcom `nlcom' (b`k':_b[/b`k'])
				}
		
				//if "`intconstraint'"=="" {
					loc b0=subinstr("`b0'","{lnshift}","_b[/lnshift]",.)
					loc b0=subinstr("`b0'","{shift}","_b[/shift]",.)
					forvalues k=1/`=`polynomial'' {
						loc b0=subinstr("`b0'","{b`k'}","_b[/b`k']",.)
					}
					forvalues b=1/`=`H'+`L'' {
						loc b0=subinstr("`b0'","{bunch`b'}","_b[/bunch`b']",.)
					}
					loc nlcom `nlcom' (b0: `b0')
					
					//}
				//else loc nlcom `nlcom' (b0: _b[/b0])
				
				forvalues k=1/`=`polynomial'' {
					if "`log'"=="" loc nlcom `nlcom' (g`k': _b[/b`k']*(1+`shift')^(`=`k'+1'))
					else {
						loc str _b[/b`k']
						if (`polynomial'>`k') {
							forvalues n=`=`k'+1'/`=`polynomial'' {
							loc str `str' +_b[/b`n']*comb(`n'+`k',`n')*ln(1+`shift')^`n'
						}
						}
					loc nlcom `nlcom' (g`k': `str')
					}
				}
				if "`log'"=="" loc nlcom `nlcom' (g0: (`b0')*(1+`shift')) 
				else {
					loc str `b0'
					forvalues n=1/`=`polynomial'' {
							loc str `str' +_b[/b`n']*ln(1+`shift')^`n'
						}
						loc nlcom `nlcom' (g0: `str')
				}
				
				if "`positiveshift'"=="positiveshift" loc nlcom `nlcom' (shift: exp(_b[/lnshift]))
				else loc nlcom `nlcom' (shift: _b[/shift])
				loc B=0
				forvalues b=1/`=`H'+`L'' {
					loc B `B' + _b[/bunch`b']
					}
				loc nlcom `nlcom' (number_bunchers: `B')
				loc m `b0'
				forvalues k=1/`polynomial' {
					loc m `m' + _b[/b`k']*`cutoff'^(`k')
				}
				loc nlcom `nlcom' (excess_mass: (`B')/(`m'))
			
				if "`log'"!="log" loc nlcom `nlcom' (marginal_response: `shift'*`cutoff')
				else loc nlcom `nlcom' (marginal_response: ln((1+`shift')*`cutoff')-ln(`cutoff'))
				
				if "`t0'"!=""&"`t1'"!="" loc nlcom `nlcom' (elasticity: ln(1+`shift')/(ln(1-`t0')-ln(1-`t1')))
				`noisily' nlcom `nlcom', post
				
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
					loc names `names' shift:lnshift
					forvalues k=1/`polynomial' {
						loc names `names' h0:``k''
					}
				}
				
				tempname b V
				mat `b'=e(b)
				mat `V'=e(V)
				
				
				
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
				estadd scalar chi2=`chi2'
				estadd scalar p_chi2=`p'
				estadd scalar df_chi2=`df_chi2'
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
					di "Likelihood ratio test of model assumptions: {col 49}Chi2(`df_chi2') test statistic {col 72}`: di %12.4f `chi2''"
					di "{col 49}p-value {col 72}`: di %12.4f `p''"
					di "{hline 83}"
				}
				
			}
				

		end
		
		
cap prog drop varcorrect
program define varcorrect, rclass
syntax anything,  fw(varname) numbins(real) [nl]
	
	qui {
		gettoken y anything: anything
		tempvar res
		tempname g V
		if "`nl'"=="nl" {
			predictnl `res'=predict() in 1/`=`numbins'', g(`g')
			replace `res'=-`res'
			mata: st_matrix("`V'",varcorrect(st_data((1,`=`numbins''),"`g'*"),st_data((1,`=`numbins''),"`fw'"),st_data((1,`=`numbins''),"`res'"),0))
			}
		else {
			predict `res' in 1/`numbins'
			replace `res'=-`res'
			mata: st_matrix("`V'",varcorrect(st_data((1,`=`numbins''),"`anything'"),st_data((1,`=`numbins''),"`fw'"),st_data((1,`=`numbins''),"`res'"),1))
		}
		return matrix V=`V'
	}
	end

mata:
	
function varcorrect(real matrix X, real matrix fw, real matrix e, real scalar addcons)
				{
				B=length(fw)
				N=sum(fw)
				if (addcons==1) X=X,J(rows(X),1,1)
				meat=0
				for (i=1; i<=B; i++) {
					e[i]=e[i]+N
					meat=meat:+fw[i]:*(X' * e * e' * X)
					e[i]=e[i]-N
				}
				return(invsym(quadcross(X,X):*N)*meat*invsym(quadcross(X,X):*N))
				}
end