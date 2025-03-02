data neph;
	length gender $5 age_g monila_g mumps_g ppd_g pha_g sksd_g $10;
	infile '~/BIOS680/data/nephroma.dat';
	
	input ptid age gender smonth sday syear
	response dmonth dday dyear status 
	monila mumps ppd pha sksd;
	
	*gender subgroup;
	gender_g = (gender='F'); 
	
	*age subgroup;
	age_g = (age>=60);
	
	*monila subgroup;
	monila_g = (monila>0);

	*mumps subgroup;
	mumps_g = (mumps>0);
	
	*ppd subgroup;
	ppd_g = (ppd>0);
	
	*pha subgroup;
	pha_g = (pha>0);
	
	*sksd subgroup;
	sksd_g = (sksd>0);
	
	* calculate survival time (days);
	startdate = mdy(smonth, sday, syear);
	enddate = mdy(dmonth, dday, dyear);
	format startdate date9.;
	format enddate date9.;
	
	dur = enddate - startdate; 

run;



* product limit (KM) estimates and log-rank test between subgroups;
* survival and hazard plots;

%macro km (var=);

ods exclude LogrankHomCov WilcoxonHomCov TaroneHomCov PetoHomCov ModPetoHomCov FlemingHomCov;
title "&var survival data and log-rank strata test";
proc lifetest data=neph method=KM plots=s(test);
	strata &var / test=all;
	time dur * status(0);
run;
title;

%mend km;

%km(var=gender_g);
%km(var=age_g);
%km(var=monila_g);
%km(var=mumps_g);
%km(var=ppd_g);
%km(var=pha_g);
%km(var=sksd_g);


* Non-parametric, log-rank strata test;
* Test association of stratified covariates with survival time
******************************************
Hypothesis test:
* H0: The subgroups within a stratitifed covariate are not different in survival time. 
The stratified covariate is not associated with survival time. 
* HA: The subgroups within a stratified covariate are different in survival time. 
The stratified covariate is associated with survival time.

Result:
None of the stratified covariates were significantly associated with surivival time at 0.05 significance level.
Only mumps skin test was nearly significant with chi-squared statistic 
3.65 df 1 (p-value=0.0561). This covariate was subgrouped into 2 levels: 
0 mm and >0 mm for the mumps skin test. 
If the covariate had been significant, 
the mumps subgroup greater than 0mm values had a longer survival time 
than mumps subgroup equal to 0mm.


Interpretation:
We fail to reject H0 for all stratified covariates at alpha=0.05 level.
Age was not a significant covariate in the log-rank stratified test,
but was significant in the non-stratitified log rank test. 
The reason may be due to loss of statistical power by stratifying age into two subgroups: 
high-age for >=60 years and low-age for <60 years;


title 'Non-parametric test on quantitative vars';
ods select CensoredSummary LogUniChiSq LogForStepSeq;
proc lifetest data=neph;
	time dur*status(0);
	test gender_g age monila mumps ppd pha sksd;
run;
title;

* Non-parametric, log-rank univariate test;
* Test association of non-stratified covariates with survival time
***************************************
Hypothesis test:
* H0: The covariate is not associated with survival time when on its own.
* HA: The covariate is associated with survival time when on its own.

Result:
Only age is significantly associated with survival time, as a continuous variable on its own.
The Log-Rank test statistic for age is -71.98, so there is a negative association between age and surivival time.

Interpretation: 
This indicates that older-aged patients have a shorter survival time than younger-aged patients.
Younger-aged patients is the only subgroup to have significant longer surivival times.
The chi-square statistic is 5.23 1 df (p-value=0.022), which is statistically significant at alpha=0.05, 
so we reject H0 and conclude that age is a significant covariate for surivival time. 



***************************************************************;

title 'Cox proportional model for 7 covariates' ;
proc phreg data=neph;
	class gender;
	model dur*status(0) = age gender monila mumps ppd pha sksd;
run;
title;



* Cox proportional model for 7 covariates;
***************************************
Global Wald / LRT

Hypothesis test: 
* H0: All covariate coefficients are equal to 0.
* HA: At least one of the covariate coefficients is not equal to 0 
or at least one covariate has a significant effect on survival time.

Result: 
* Fitting all 7 covariates in the model results in a chi-squared statistic above 10 df 7 (p-value=0.18).
We fail to reject H0 and conclude that fitting 7 covariates does not show significant effects on survival time. 

Interpretation: 
* Overfitting too many covariates does not show which covariates are most significant. 
We refit the model selecting fewer covariates (age and mumps) based on the non-parametric results.; 


title 'Cox proportional model for 2 covariates' ;
proc phreg data=neph;
	model dur*status(0) = age mumps;
run;
title;

* Global Wald / LRT

Hypothesis test:  
* H0: All covariate coefficients are equal to 0.
* HA: At least one of the covariate coefficients is not equal to 0 
or at least one covariate has a significant effect on survival time.

Result: 
* Fitting age and mumps in the Cox proportional model results in LRT chi-squared statistic 
at 6.96 df 2 (p-value=0.031) and Wald chi-squared statistic 7.22 df 2 (p-value=0.027). 
We reject H0 and conclude that at least one of the covariate coefficients is not equal to 0. 


Interpretation: 
* Next, we test each covariate, age and mumps using Cox proportional hazard test.

*******************************************
* Cox proportional hazard test: 

Hypothesis test:
* H0: The covariate is not associated with survival time after adjusting for the other covariate in the model.
* HA: The covariate is associated with survival time after adjusting for the other covariate in the model.

Result: 
Only age is significantly associated with survival time after adjusting for mumps skin test.
The chi-squared test statistic for age is 4.57 df 1 (p-value=0.0324).
We reject H0 and conclude that age has a significant effect on surivival time, after adjusting for mumps skin test.

Age was also a significant covariate in the log-rank univariate test.

Interpretation: 
The parameter estimate for age is 0.06785, and the hazard ratio is exp(0.06785) = 1.070. 
For every 1-year increase in age of patient, the hazard for death goes up by an estimated 7.0 percent.




