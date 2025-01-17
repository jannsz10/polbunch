{smcl}
{cmd:help polbunchgendata}
{hline}

{title:Title}

{p2colset 5 14 16 2}{...}
{p2col:{cmd:polbunchgendata} {hline 2}}Simulates data from the iso-elastic labor supply model in a setting with a kink.{cmd: polbunch}
{p_end}
{p2colreset}{...}


{marker syntax}{...}
{title:Syntax}

{p 10 15 2}
{cmd:polbunchgendata} newvarname {cmd:,} obs(#) cutoff(#) el(#) t0(#) t1(#) [distribution(string) log]

{pstd}
{cmd:polbunchgendata} generates data from the iso-elastic labor supply model with a tax kink at {cmd: cutoff} and a counterfactual distribution 
of earnings specified in {opt: distribution}.

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{synopt:{opt obs(#)}} Number of observations.{p_end}
{synopt:{opt t0(#)}} Tax rate below the kink.{p_end}
{synopt:{opt t1(#)}} Tax rate above the kink.{p_end}
{synopt:{opt cutoff(#)}} Kink point.{p_end}
{synopt:{opt distribution(string)}} A string that will generate a valid random distribution of earnings, for instance rbeta(2,5). Also allowed are triangular(a,b,c), with triangular(0,3,0) the default.{p_end}
{synopt:{opt log}} Specifies that the earnings variable is in logs.T{p_end}
{synoptline}

{marker Author}{...}
{title:Author}

{pstd}Martin Eckhoff Andresen{p_end}
{pstd}University of Oslo{p_end}
{pstd}Department of Economics{p_end}
{pstd}Oslo, Norway{p_end}
{pstd}martin.eckhoff.andresen@gmail.com{p_end}

{marker also_see}{...}
{title:Also see}

{p 7 14 2}
{helpb polbunch} for help on the main command {cmd: polbunch}.{p_end}
