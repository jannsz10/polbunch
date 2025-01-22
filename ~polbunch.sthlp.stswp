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
{synopt:{opt noi:sily}} display intermediate output such as iteration logs.{p_end}
{synopt:{opt nodots}} suppress dots indicating bootstrap progress{p_end}
{synopt:{opt notest}} do not test model assumptions{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:polbunch} performs polynomial bunching estimation using OLS and nonlinear least squares.

There are two different ways to use polbunch. If one variable is specified, polbunch interprets this as data containing individual level data on earnings. 
The option bw() must then be specified, and polbunch collapses the data to bins before performing bunching estimation. During this step, the leftmost and 
rightmost bins are dropped if polbunch finds that there is limited coverage in these bins, unless the option nodrop is specified. If two variables are 
specified, polbunch interprets this as data being already binned, and take the first variable to contain the bin counts and the second to containt the 
bin midpoints. Option bw() cannot be specified, and will be determined by polbunch from the bin widths. 

Polbunch implements four different bunching estimators, but beware that these are supported for comparison to alternative estimators in the litterature, 
and besides the default, the other alternatives are neither consistent with theory, unbiased or efficient. The default choise of estimator(3) implements 
the theoretically consistent and efficient polynomial bunching estimator of Andresen (2025). Like most bunching estimators, this relies on very strong 
structural assumptions, including the correct specification of the iso-elastic labor suppply model and a counterfactual density that is well approximated 
by the specified polynomial.

Alternative estimators include 
- 0 for the unbiased, but inefficient estimator that only relies on data to the left of the cutoff to estimate the counterfactual density.
- 1 for the biased and inconsistent estimator that ignores the intensive margin response to the higher tax rate above non-bunching agents above the cutoff
- 2 for the adjustment procedure of Chetty et. al. (2011).
- 3 for the default, the efficient and theoretically consistent estimator of Andresen (2025)

After having constructed the binned data, polbunch proceeds by estimating a regression of the bin counts on separate polynomials below and above the cutoff
 and a series of fixed effects for each bin in the excluded region. If estimator(0) is specified, the results are posted based on this regression.

If one of the more restrictive estimators 1-3 are specified, polbunch then tests the restrictions imposed on this unrestricted model by that estimator, 
and reports the Chi2 test statistic and associated p-value, unless the option "notest" is specified. This may be interpreted as a joint test of whether 
the restrictions implicitly imposed by the estimator is consistent with the coefficients on the separately estimated polynomials from the unrestricted model, 
and the parametric polynomial assumption itself. Polbunch then estimates the corresponding restricted model using nonlinear least squares, based on initial
 values from the linear regression.

If the option notransform is specified, polbunch reports the parameters of the estimating equation, which depends on the estimator. Otherwise, and as a 
default, polbunch transforms the reported estimates and report the parameters of h0 (the density of earnings under the low tax), h1 (the density of 
earnings under the high tax) and various bunching parameters. These include the number of bunchers, the excess mass, the proportional shift parameter, 
the response of the marginal buncher and, if t0() and t1() are specified, the elasticity.

Inference for the reported parameters is based on the analytical standard errors (the default, if bootreps(1)) or the binned bootstrap procedure 
(if bootreps>=2) developed by Andresen (2025), unless the bootreps(0) option is specified, which omits the calculation of variance of the estimated 
results. This may be useful if the user wish to perform a regular bootstrap procedure with the {cmd: bootstrap} prefix or during simulation exercises
 where only the parameters are of interest. 

See {cmd: polbunchplot} for a companion command to plot estimated densities and results from polbunch.

{marker options}{...}
{title:Details on selected options not discussed above}

{phang}
{cmd: constant} specifies that polbunch should use the constant approximation to estimate the response of the marginal buncher
 and the elasticity. This is innocuous if the true elasticity is small, the tax change is marginal or the true density is flat, 
 but otherwise may induce bias. The default is to solve the integral under the estimated density for the response of the marginal 
 buncher that exactly fits the number of bunchers under the counterfactual, which do not suffer from this approximation bias. 
 With this default method, however, analytical standard errors cannot be derived for the response of the margina buncher or the 
 elasticity, so that polbunch do not report standard errors for these parameters and the user needs to rely on the bootstrap 
 using the bootreps(#) option with #>2.

{phang}
{cmd: nopositiveshift} specifies that there should be no restriction on the shift parameter during estimation. By default, 
this is restricted to be nonnegative, implying that the elasticity is also nonnegative. 


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

{phang2}{cmd:. polbunchgendata, obs(10000) t0(0.2) t1(0.6) el(0.4) cutoff(1)}{p_end}

{pstd}
Estimate bunching using the default estimator and a correctly specified polynomial of degree 1 and a bandwidth of 0.01.{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) pol(1)}{p_end}

{pstd}
Instead use the (biased) esitmator of Chetty et. all.{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) pol(1) estimator(2)}{p_end}

{pstd}
Or the naive, biased esitmator that ignores intensive margin response above the threshold.{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) pol(1) estimator(1)}{p_end}

{pstd}
Use the binned bootstrap for inference.{p_end}
{phang2}{cmd:. polbunch z, cutoff(1) bw(0.01) pol(1) bootreps(200)}{p_end}


{marker saved_results}{...}
{title:Stored results}

{pstd}
{cmd:polbunch} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations used (after trimming unless nodrop specified){p_end}
{synopt:{cmd:e(F_mod)}}F-test of model specification, if estimator>0{p_end}
{synopt:{cmd:e(p_mod)}}p-value of test for model specification, if estimator>0{p_end}
{synopt:{cmd:e(F_df1)}}degrees of freedom, F-test of model specification{p_end}
{synopt:{cmd:e(F_df2)}}degrees of freedom, F-test of model specification{p_end}
{synopt:{cmd:e(polynomial)}}degree of polynomial used{p_end}
{synopt:{cmd:e(cutoff)}}kink position{p_end}
{synopt:{cmd:e(bw)}}bandwidth{p_end}
{synopt:{cmd:e(lower_limit)}}lower limit of excluded region{p_end}
{synopt:{cmd:e(upper_limit)}}upper limit of excluded region{p_end}
{synopt:{cmd:e(aic)}}Akaike Information Criterion{p_end}
{synopt:{cmd:e(log)}}Whether the earnings variable was specified in logs{p_end}
{synopt:{cmd:e(p_X)}}p-value for a test of observed heterogeneity -- the b_1-b_0 coefficients{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(binname)}}name of the earnings variable{p_end}
{synopt:{cmd:e(cmd)}}{cmd:polbunch}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(depvar)}}freq (if no pre-binned data) or name of the dependent variable{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}; {cmd:V}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector, including counterfactual densities and bunching parameters{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of {cmd:b}.{p_end}
{synopt:{cmd:e(table)}}table of frequencies.{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{marker references}{...}
{title:References}

{phang}
Andresen, Martin E. (2025). "A better polynomial bunching estimator", working paper.

{phang}
Kleven, Henrik Jacobsen (2016). "Bunching", Annual Review of Economics

{phang}
Saez, Emmanuel (2010), "Do Taxpayers bunch at kink points?", American Economic Review

{phang}
Chetty, Raj, Friedman, John N., Olsen, Tore and Luigi Pistaferri (2011), "Adjustment Costs, Firm Responses, and Micro vs. Macro Labor Supply Elasticities: Evidence from Danish Tax Records", Quarterly Journal of Economics


{title:Thanks for citing pol as follows}

{pstd}
Andresen, Martin E., (2025). "POLBUNCH: Stata module for the polynomial bunching estimator.". This version VERSION_DATE.{p_end}

{pstd}
where you can check your version date as follows:{p_end}

{phang2}{cmd:. which polbunch}{p_end}

	 
{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}University of Oslo{p_end}
{pstd}Department of Economics{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{marker also_see}{...}
{title:Also see}

{p 4 14 2}
Development version: net install mtefe, from("https://raw.githubusercontent.com/martin-andresen/polbunch/master"){p_end}

{p 7 14 2}
Help:  {helpb polbunchplot}, {helpb polbunchgendata}{p_end}
