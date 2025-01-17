{smcl}
{cmd:help polbunch}{right: ()}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:polbunch} {hline 2}}Efficient and theoretically consistent polynomial bunching estimation
{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 8 13 2}
{cmd:polbunch} [{freqvar}] {zvar} {ifin}{cmd:,} cutoff(#) [ {it:options}]

			NOIsily ///
			nodots /// suppress dots for bootstrap progress
			]
			
{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Options}
{synopt:{opt cut:off(#)}} is required and specifies the kink point. {p_end}
{synopt:{opt pol:ynomial(#)}} specifies the degree of the polynomial, default is 7. {p_end}
{synopt:{opt est:imator(#)}} specifies the estimator - the default is 3, the efficient and consistent polynomial bunching estimator. {p_end}
{synopt:{opt boot:reps(#)}} specifies the number of bootstrap repetitions. The default of 1 implies analytic standard errors, 0 implies no standard errors.{p_end}
{synopt:{opt log}} specifies that the earnings variable is in logs{p_end}
{synopt:{opt t0}} specifies the linear tax rate below the cutoff{p_end}
{synopt:{opt t1}} specifies the linear tax rate above the cutoff{p_end}
{synopt:{opt lim:its(numlist)}} specifies how many bins are excluded below the cutoff (1st number) and how many above (2nd number), the default is 1 0.{p_end}
{synopt:{opt constant}} uses the constant approximation to calculate the response of the marginal buncher and the elasticity.{p_end}
{synopt:{opt notransform}}reports the coefficients of the estimation model rather than transforming these to bunching estimates.{p_end}
{synopt:{opt nopositiveshift}}does not restrict the shift parameter to be positive.{p_end}
{synopt:{opt initvals(string)}}specifies initial values for the NLS procedure. By default polbunch finds reasonable feasible values from an initial OLS regression.{p_end}
{synopt:{opt nodrop}} implies that polbunch does not drop the leftmost and rightmost bin when these appear to be cut by sample selection.{p_end}
{synopt:{opt noi:sily}} display intermediate regression output such as iteration logs.{p_end}
{synopt:{opt nodots}} suppress dots indicating bootstrap progress

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunch} performs polynomial bunching estimation using OLS and nonlinear least squares.

There are two different ways to use polbunch. If one variable is specified, polbunch interprets this as data containing individual level data on earnings. The option bw() must then be specified, and polbunch collapses the data to bins before performing bunching estimation. During this step, the leftmost and rightmost bins are dropped if polbunch finds that there is limited coverage in these bins, unless the option nodrop is specified. If two variables are specified, polbunch interprets this as data being already binned, and take the first variable to contain the bin counts and the second to containt the bin midpoints. Option bw() cannot be specified, and will be determined by polbunch from the bin widths. 

Polbunch implements four different bunching estimators, but beware that these are supported for comparison to alternative estimators in the litterature, and besides the default, the other alternatives are neither consistent with theory, unbiased or consistent. The default choise of estimator(3) implements the theoretically consistent and efficient polynomial bunching estimator of Andresen (2025). Like most bunching estimators, this relies on very strong structural assumptions, including the correct specification of the iso-elastic labor suppply model and a counterfactual density that is well approximated by the specified polynomial.

Alternative estimators include 
- 0 for the unbiased, but inefficient estimator that only relies on data to the left of the cutoff to estimate the counterfactual density.
- 1 for the biased and inconsistent estimator that ignores the intensive margin response to the higher tax rate above non-bunching agents above the cutoff
- 2 for the adjustment procedure of Chetty et. al. (2011).
- 3 for the default, the efficient and theoretically consistent estimator of Andresen (2025)

After having constructed the binned data, polbunch proceeds by estimating a regression of the bin counts on separate polynomials below and above the cutoff and a series of fixed effects for each bin in the excluded region. If estimator(0) is specified, the results are posted based on this regression.

If one of the more restrictive estimators 1-3 are specified, polbunch estimates the corresponding model using nonlinear least squares, based on initial values from the linear regression.

If the option notransform is specified, polbunch reports the parameters of the estimating equation, which depends on the estimator. Otherwise, and as a default, polbunch transforms the reported estimates and report the parameters of h0 (the density of earnings under the low tax), h1 (the density of earnings under the high tax) and various bunching parameters. These include the number of bunchers, the excess mass, the proportional shift parameter, the response of the marginal buncher and, if t0() and t1() are specified, the elasticity.

Inference for the reported parameters is based on the analytical standard errors (the default, if bootreps(1)) or the binned bootstrap procedure (if bootreps>=2) developed by Andresen (2025), unless the bootreps(0) option is specified, which omits the calculation of variance of the estimated results. This may be useful if the user wish to perform a regular bootstrap procedure with the {cmd: bootstrap} prefix or during simulation exercises where only the parameters are of interest. 

See {cmd: polbunchplot} for a companion command to plot estimated densities and results from polbunch.

{marker options}{...}
{title:Details on options}

{phang}
{opt polynomial(#)} specifies the k(p) or k_j(p) functions to be polynomials
of degree {it:#}. Required for the parametric and semiparametric polynomial
models.  Note that {cmd:mtefe} specifies the degree of the MTE directly,
rather than the degree of K(p) like {cmd:margte}, so a polynomial of L+1 in
{cmd:margte} will be equivalent to L in {cmd:mtefe}.


{marker examples}{...}
{title:Examples}

{pstd}
To provide examples, the {cmd:polbunchgendata} program simulates data on earnings
based on a user-specified distribution of earnings in the low tax-regime and the 
iso-elastic labor supply model.

{pstd}
Generate simulated data of 10,000 observations of log wages that would be distributed
according to a triangular(0,3,0) distribution in the low-tax regime, where the
earnings follow the optiomal choiced from the iso-elastic labor supply model, the 
elasticity is 0.4 and the tax rates below and above the cutoff is t0=0.2 and t1=0.6{p_end}

{phang2}{cmd:. polbunchgendatagendata, obs(10000) t0(0.2) t1(0.6) el(0.4) cutoff(1)}{p_end}

{pstd}
The parametric normal model, fit using local IV{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol)}{p_end}

{pstd}
The parametric polynomial of degree 1, estimated using the separate
approach{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), separate pol(1)}{p_end}

{pstd}
The polynomial model of degree 2, using a linear probability model for the
propensity score{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), polynomial(2) link(lpm)}{p_end}

{pstd}
The parametric normal model, fit using maximum likelihood{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), mlikelihood}{p_end}

{pstd}
The semiparametric polynomial model with 50 bootstrap replications {p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), pol(1) bootreps(50)}{p_end}

{pstd}
The joint normal model where the fixed effects are restricted to be the same in
the treated and untreated state (note: violation of exclusion) {p_end}
{phang2}{cmd:. mtefe lwage exp exp2 (col=distCol), restricted(i.district)}{p_end}

{pstd}
The semiparametric model, using a smaller evaluation grid to save computational time{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), semiparametric gridpoints(100)}{p_end}

{pstd}
Calculating the PRTE for a policy that mandates a maximum distance to college
of 40 miles, normal model{p_end}
{phang2}{cmd:. probit col distCol exp exp2 i.district}{p_end}
{phang2}{cmd:. rename distCol tempdistcol}{p_end}
{phang2}{cmd:. gen distCol=min(40,tempdistcol)}{p_end}
{phang2}{cmd:. predict double p_prte}{p_end}
{phang2}{cmd:. drop distCol}{p_end}
{phang2}{cmd:. rename tempdistcol distCol}{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), prte(p_prte)}{p_end}
{phang2}{cmd:. mtefeplot, prte}{p_end}

{pstd}
The polynomial model of degree 2 with two splines at 0.25 and 0.75, fit using
the local IV{p_end}
{phang2}{cmd:. mtefe lwage exp exp2 i.district (col=distCol), polynomial(2) splines(0.25 0.75)}{p_end}


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:mtefe} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(IV)}}regular two-stage least-squares estimate{p_end}
{synopt:{cmd:e(polynomial)}}degree of polynomial in polynomial models{p_end}
{synopt:{cmd:e(p_U)}}p-value for a test of unobserved heterogeneity -- the coefficients on k(p){p_end}
{synopt:{cmd:e(p_X)}}p-value for a test of observed heterogeneity -- the b_1-b_0 coefficients{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:mtefe}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(title2)}}secondary title in estimation output{p_end}
{synopt:{cmd:e(method)}}estimation method{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}; {cmd:V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector, including MTE and treatment parameters{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:b}, empty if using semiparametric methods and no bootstrap{p_end}
{synopt:{cmd:e(trimminglimits)}}trimming limits if using the trimsupport() option. {p_end}
{synopt:{cmd:e(support)}}vector of common support over p{p_end}
{synopt:{cmd:e(mte)}}vector of MTE estimates{p_end}
{synopt:{cmd:e(mtexs_ate)}}vector of x values for which the main MTE and ATE
were calculated{p_end}
{synopt:{cmd:e(mteatt)}}MTE curve for the treated individuals{p_end}
{synopt:{cmd:e(weightsatt)}}estimated weights at each point of support for the ATT{p_end}
{synopt:{cmd:e(mtexs_att)}}vector of x values for which the ATT was calculated{p_end}
{synopt:{cmd:e(mteatut)}}MTE curve for the untreated individuals{p_end}
{synopt:{cmd:e(weightsatut)}}estimated weights at each point of support for the ATUT{p_end}
{synopt:{cmd:e(mtexs_atut)}}vector of x values for which the ATUT was calculated{p_end}
{synopt:{cmd:e(mtelate)}}MTE curve for the compliers{p_end}
{synopt:{cmd:e(weightslate)}}estimated weights at each point of support for the LATE{p_end}
{synopt:{cmd:e(mtexs_atut)}}vector of x values for which the LATE was calculated{p_end}
{synopt:{cmd:e(mteprte)}}MTE curve for the policy compliers if specified in {cmd:prte()}{p_end}
{synopt:{cmd:e(weightsprte)}}estimated weights at each point of support for the PRTE (if specified in {cmd:prte()}){p_end}
{synopt:{cmd:e(mtexs_atut)}}vector of x values for which the PRTE was calculated (if specified in {cmd:prte()}){p_end}
{synopt:{cmd:e(dkdp)}}estimates of k(u) at the points of support{p_end}
{synopt:{cmd:e(weightsmprte1)}}estimated weights at each point of support for the MPRTE1 parameter{p_end}
{synopt:{cmd:e(weightsmprte2)}}estimated weights at each point of support for the MPRTE2 parameter{p_end}
{synopt:{cmd:e(weightsmprte3)}}estimated weights at each point of support for the MPRTE3 parameter{p_end}
{synopt:{cmd:e(tescales)}}scales of the non-marginal treatment effect parameters in situations with limited support{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{pstd}
If {cmd:semiparametric} is specified, the following matrices are stored:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(bandwidth)}}bandwidths used in local polynomial smooths, alternatively as {cmd:e(bandwidth1)} and {cmd:e(bandwidth0)} if the separate approach is used{p_end}
{synopt:{cmd:e(degree)}}degree in local polynomial smooth{p_end}

{pstd}
If {cmd:separate} or {cmd:mlikelihood} is specified, the following matrices
are stored:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(Y1)}}average outcomes in the treated state{p_end}
{synopt:{cmd:e(Y0)}}average outcomes in the untreated state{p_end}

{pstd}
Note: {cmd:bootstrap} stores its own output in {cmd:e()} as well.  See
{helpb bootstrap##saved_results:bootstrap}.

{pstd}
Note: Various options store results as variables as well. These include
{cmd:saveweights()}, where the variable prefix is user specified;
{cmd:savepropensity}, where the new variable name is user specified; and
{cmd:savekp}, where the variable names are predetermined to work with
postestimation predictions.




{marker references}{...}
{title:References}

{phang}
Brave, S., and T. Walstrum. 2014. {browse "http://www.stata-journal.com/sjpdf.html?articlenum=st0331":Estimating marginal treatment effects using parametric and semiparametric methods}. {it:Stata Journal} 14: 191-217.

{phang}
Brinch C., M. Mogstad, and M. Wiswall. 2017. Beyond LATE with a discrete instrument. {it:Journal of Political Economy} 125: 985-1039.

{phang}
Heckman, J. J., and E. J. Vytlacil. 2007a. Econometric evaluation of social
programs, part I: Causal models, structural models and econometric policy
evaluation. In {it:Handbook of Econometrics}, vol. 6B, ed. J. J. Heckman and E.
E. Leamer, 4779-4874. Amsterdam: Elsevier.

{phang}
------. 2007b. Econometric evaluation of social programs, part II: Using the
marginal treatment effect to organize alternative econometric estimators to
evaluate social programs, and to forecast their effects in new environments. In
{it:Handbook of Econometrics}, vol. 6B, ed. J. J. Heckman and E. E. Leamer,
4875-5143. Amsterdam: Elsevier.

{phang}
Lokshin, M. and Z. Sajaia. 2004.
{browse "http://www.stata-journal.com/sjpdf.html?articlenum=st0071":Maximum likelihood estimation of endogenous switching regression models}.
{it:Stata Journal} 4: 282-289.


{title:Thanks for citing mtefe as follows}

{pstd}
Andresen, Martin E., (2018). "MTEFE: Stata module to estimate marginal treatment effects". This version VERSION_DATE.{p_end}

{pstd}
where you can check your version date as follows:{p_end}

{phang2}{cmd:. which mtefe}{p_end}

	 
{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}Statistics Norway{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{pstd}
This program is inspired by {helpb margte} by Brave and Walstrum (2014), who
deserve a huge thanks.  The maximum likelihood evaluator is adapted from
the {helpb movestay} command by Lokshin and Sajaia. Edwin Leuven, Katrine Løken, 
Magne Mogstad also deserve thanks for various help and bug reports.


{marker also_see}{...}
{title:Also see}

{p 4 14 2}
Development version: net install mtefe, from("https://raw.githubusercontent.com/martin-andresen/mtefe/master"){p_end}

{p 4 14 2}
Article:  {it:Stata Journal}, volume 18, number 1: {browse "http://www.stata-journal.com/article.html?article=st0516":st0516}{p_end}

{p 7 14 2}
Help:  {helpb mtefeplot}, {helpb locpoly3}, {helpb nearmrg} (required for {cmd:gridpoints()}), {helpb movestay}, {helpb margte} (if installed){p_end}
