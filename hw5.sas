
data burn;
	length sev $4;
	infile '~/BIOS680/data/burn.txt';
	input obs z1 z2 z3 z4 z5 z6 z7 z8 z9 z10 z11 
	t1 delta1 t2 delta2 t3 delta3;
		
		* separate burn severity by grade level (categorical);
		if z4>= 0 and z4 < 15 then sev = 'low';
		else if z4 >=15 and z4 < 70 then sev = 'med';
		else sev = 'high';
run;


* 2.
The burn.dat has data on the patients with burn-wound infections. 
Patients were followed until a staphylococcus infection occurred. 
Of primary interest in the study is a comparison
of two methods of body cleansing: 
a routine bathing care method and a body cleansing
method using 4% chlorhexidine gluconate. 
Create a burn severity variable based on the Percentage of total surface area burned 
(0-14.99%, 15%-69.99% and 70%+). 
In the following tests, clearly state your H0 and HA. and using alpha = 0:05.;

* (2a) Plot the K-M curves for the two treatments and test if the two treatments are the
same in the overall sample.

Log-Rank Test for treatment groups 1=routine 2=special treatment
p-value = 0.0515 , fail to reject H0

H0: Both routine bathing care and 4% chlorhexidine gluconate treatment have equal staph infection probability. 
	  S(t1) = S(t2)
Ha: Routine bathing care differs in staph infection probability from 4% chlorhexidine gluconate treatment cleansing.
  	  S(t1) ^= S(t2) or S(t1) = c*S(t2) for constant c^=1;


title '2a) Plot K-M curves of two treatments';
proc lifetest data=burn;
	ods select HomTests SurvivalPlot;
	strata z1;
	time t3 * delta3(0);
run;
title;

* (2b) Generate K-M curves for the two treatments for each of the 3 strata defined by the
burn severity. Plot them on the same page. ;


title '2b) Plot K-M curves of two treatments for each of 3 strata (burn severity)';
proc lifetest data=burn;
	ods select HomTests SurvivalPlot;
	strata sev/ group=z1;
	time t3 * delta3(0);
run;
title;

*Do a stratified test on the treatments
based on the burn severity.
Stratified log-rank test 
p-value = 0.0719 , fail to reject H0

H0: h(t1|severity) = h2(t2|severity)
Ha: h(t1|severity) ^= h2(t2|severity) or h(t1|severity) = c*h2(t2|severity) for constant c^=1;

* (2c) Comment on a and b.
Adjusting for burn severity, the two treatments do not differ in staph infection hazard rate, 
and thus do not differ in staph infection probability. 
Also without adjustment, the two treatments do not differ in staph infection probability.
 
* (2d) Test whether the hazards increase with the burn severity.
Test for trend, burn severity levels 1=low, 2=med, 3=high

H0: h1(t) = h2(t) = h3(t) 
HA: h1(t) <= h2(t) <= h3(t) 

p-value = 0.3956, fail to reject H0.
There is not sufficient evidence that hazard of staph infection increases with increasing burn severity level;

proc format;
    value $fmt
        'low' = '1L'
        'med' = '2M'
        'high' = '3H';
run;

title '2d) Test of trend';
proc lifetest data=burn;
	ods select HomTests SurvivalPlot;
	strata sev / order=formatted;
	time t3 * delta3(0);
	format sev $fmt.;
run;
title;
  
  
  
***************************

* 3.;
data oc;
	infile '~/BIOS680/data/oc.txt';
	input age grade stage disease trmt time event;
	if (grade^=9); * 9 = missing; *remove missing grade values;
	g2 = (grade=2);
	g3 = (grade=3);
	g_1 = grade -1;
	p32g_1 = p32*g_1;
run;

proc freq data=oc;
	table grade / missing list;
run;

*(3a) Estimate the survival function for each grade. Make plots to check graphically
whether the exponential or the Weibull distributions seem reasonable. Comment.;




title ' 3a';
proc lifetest data=oc plots=(hazard, s) outs=outkm1 alpha=0.05;
	ods select SurvivalPlot HazardPlot;
	strata grade;
	time time * event(0);
	where grade = 1;
run;
title;

title ' 3a';
proc lifetest data=oc plots=(hazard, s) outs=outkm2 alpha=0.05;
	ods select  SurvivalPlot HazardPlot;
	strata grade;
	time time * event(0);
	where grade = 2;
run;
title;

title ' 3a';
proc lifetest data=oc plots=(hazard, s) outs=outkm3 alpha=0.05;
	ods select  SurvivalPlot HazardPlot;
	strata grade;
	time time * event(0);
	where grade = 3;
run;
title;

title ' 3a Survival function estimates (KM) for grade1';
proc print data=outkm1;
run;
title;

title ' 3a Survival function estimates (KM) for grade2';
proc print data=outkm2;
run;
title;

title ' 3a Survival function estimates (KM) for grade3';
proc print data=outkm3;
run;
title;
*/

*
grade 1 = no distribution, 1 extreme hazard value due to a failure at a single timepoint
grade 2 = Weibull lambda > 0, gamma < 0, monotonically dec
grade 3 = non-monotonic

best prognosis of recurrence: grade 1, low hazard, recurrence rate is low
worst prognosis of recurrence: grade 3, high hazard, recurrence rate is high;

* (3b)
Fit an exponential model separately within each grade. Report the estimated hazards
with their standard errors and 95% condence interval. Report log-hazards with their
standard errors and the 95% condence intervals. Interpret the hazard estimates.;


title '3b class grade';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	class grade;
	model time*event(0) = grade 
	/ distribution=exponential;
	probplot;
run;
title;

title '3b grade1';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	where grade=1;
	model time*event(0) = grade 
	/ distribution=exponential;
	probplot;
run;
title;

title '3b grade2';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	where grade=2;
	model time*event(0) = grade 
	/ distribution=exponential;
	probplot;
run;
title;

title '3b grade3';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	where grade=3;
	model time*event(0) = grade 
	/ distribution=exponential;
	probplot;
run;
title;

* (3c)
Under the exponential model, test the null hypothesis that the three grades have
the same distribution of time to recurrence. Set the type I error at 0.05. Report
p-values from i) the likelihood ratio test and ii) the Wald test (use SAS with class
statement).;

* LRT test statistic:
/* Full models with grades, -2LL = 256.852*/
ods output FitStatistics=full;
title '3c';
proc lifereg data=oc;
	ods select FitStatistics;
	class grade;
    model time*event(0) = grade / dist=exponential;
run;

/* Reduced model without grades, -2LL = 281.683 */

proc lifereg data=oc;
	ods select FitStatistics;
    model time*event(0) = / dist=exponential;
run;

* Wald test statistic;

proc lifereg data=oc;
	ods select Type3Analysis ParameterEstimates;
	class grade;
    model time*event(0) = grade / dist=exponential;
run;
title;

* (3d);
* Fit Weibull model by each grade; 
* find estimated mean and median survival times;
* The PPlot shows that the weibull distribution does not fit 
as well as the exponential distribution (points are not linear in the plot);


title '3d class grade';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	class grade;
	model time*event(0) = grade 
	/ dist=weibull;
	probplot;
run;
title;

title '3d grade1';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	where grade=1;
	model time*event(0) = grade 
	/ dist=weibull;
	probplot;
run;
title;

title '3d grade2';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	where grade=2;
	model time*event(0) = grade 
	/ dist=weibull;
	probplot;
run;
title;

title '3d grade3';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	where grade=3;
	model time*event(0) = grade 
	/ dist=weibull;
	probplot;
run;
title;


/*title '3d';
proc lifereg data=oc;
	ods select ParameterEstimates ProbPlot;
	class grade;
	model time*event(0) = grade / distribution=weibull;
	probplot;
	output out=med predicted=median;
run;

proc sort data=med nodupkey;
	by grade;
proc print data=med;
	by grade;
run;
title;
*/


* (3e);
* For grade = 2, test the null hypothesis that the recurrence times have an exponential
distribution against the alternative that they have a Weibull distribution using the
likelihood ratio test. Set the Type I error at 0.05. Report the p-value.;
data grade2;
	set oc;
	where grade = 2;
run;

* Fit exponential dist, -2LL = 281.683;
title '3e';
proc lifereg data=oc;
	ods select ModelInfo FitStatistics;
    model time*event(0) = / dist=exponential;
run;

* Fit Weibull dist, -2LL = 264.143;
proc lifereg data=oc;
	ods select ModelInfo FitStatistics;
    model time*event(0) = / dist=weibull;
run;
title;

* (3f)
For an exponential model with hazard , for grade = 2, test the hypothesis that the
true median survival time is one year against the alternative that it is more than
one year.;
proc lifereg data=oc;
	model time*event(0) = grade / dist=exponential;
	output out=survtime predicted=median;
	where grade=2;
run;

* Wilcoxon signed rank test;
title '3f Wilcoxon signed rank test';
proc univariate data=survtime(where=(grade=2)) loccount mu0=365; /*H0: median=1*/
	ods exclude Moments BasicMeasures;
  	var median;
run;
title;
ods graphics off;

* (3g)
Take the control group as the baseline. Assume that their survival times have a
Weibull distribution. Fit a model in which treatment accelerates (or decelerates)
the time to recurrence. Interpret your parameter estimates (for example how would
you explain the results to a physician). Obtain a 95% condence interval for the
acceleration factor.;

* (3h)
Take the control group as the baseline. Assume that their survival times have a
Weibull distribution. Fit a model in which treatment scales the hazard by a factor
(the hazard ratio). Interpret your parameter estimates. Obtain a 95% condence
interval for the hazard ratio.;
title '3g-h';

proc lifereg data=oc;
	ods select ParameterEstimates ;
	model time*event(0) = trmt / dist = weibull;
run;
title;

* (4a)
Use two groups for stage: stage
3 (including 3a, 3b, 3ab, 3A, 3B, 3AB) and stage 4 (including 4a, 4b, 4A, 4B). Take stage
4 as the reference group. Fit an accelerated life model with Weibull distribution.
(4a-d) 
(a) Based on the fitted model, calculate, by hand and by SAS, the median survival time
and its 95% confidence interval for stage 4 patients.
(b) Based on the fitted model, calculate, by hand and by SAS, the median survival time
and its 95% confidence interval for stage 3 patients.
(c) Based on the fitted model, calculate the probability of surviving to 15 months for
the stage 3 patients.
(d) Based on the fitted model, calculate the probability of surviving to 15 months for
the stage 4 patients.;

data mel;
	infile '~/BIOS680/data/melanoma.txt';
	input ptid age sex fam rtime rstatus 
	stime sstatus stage $ monilia mumps ppd
	pha sksd trico;
	if upcase(stage) in ('3A', '3B', '3AB') then stages = 3;
	else stages = 4;
run;


/*
proc freq data=mel;
	table stage*stages / list missing;
run; */


proc format;
	value st
	3 = '1-stage3'
	4 = '2-stage4'
;


title '4a-b';
proc lifetest data=mel plots=(h, s);
	strata stages / order=formatted;
	time stime*sstatus(0);
run;
title;

ods trace on;
proc lifereg data=mel outest=parm covout order=formatted;
	class stages;
	model stime*sstatus(0) = stages / dist=weibull;
	output out=sfit q=0.5 predicted=median std=se;
	probplot;
	format stages st.;
run;

title '4a-b Median survival time for stage 3 or 4 patients';
proc sort data=sfit nodupkey;
	by stages;
proc print data=sfit;
	var stages median se;
	by stages;
run;
title;
  
  
  
