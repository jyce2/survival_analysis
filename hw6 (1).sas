


* 3 Fit models by each grade;
data oc;
	infile '~/BIOS680/data/oc.txt';
	input age grade stage disease trmt time event;
	if (grade^=9); * 9 = missing; *remove missing grade values;
	g2 = (grade=2); *g2 indicator = 1;
	g3 = (grade=3); *g3 indicator = 1;
	g_1 = grade -1; *g_1 continuous var as 0,1,2;
	p32g_1 = trmt*g_1; *interaction of trmt*grade;
run;

************************************************;

title '3b grade1';
proc lifereg data=oc;
	ods select ModelInfo ParameterEstimates ProbPlot;
	where grade=1;
	model time*event(0) = grade / dist=exponential;
	probplot;
	output out=csres3b1 cresidual=cri;
run;

proc lifetest data=csres3b1 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;

title;


title '3b grade2';
proc lifereg data=oc;
	ods select ModelInfo ParameterEstimates ProbPlot;
	where grade=2;
	model time*event(0) = grade / dist=exponential;
	probplot;
	output out=csres3b2 cresidual=cri;
run;

proc lifetest data=csres3b2 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;

title;

title '3b grade3';
proc lifereg data=oc;
	ods select ModelInfo ParameterEstimates ProbPlot;
	where grade=3;
	model time*event(0) = grade / dist=exponential;
	probplot;
	output out=csres3b3 cresidual=cri;
run;

proc lifetest data=csres3b3 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;

title;


title '3d grade1';
proc lifereg data=oc;
	ods select ModelInfo ParameterEstimates ProbPlot;
	where grade=1;
	model time*event(0) = grade / dist=weibull;
	probplot;
	output out=csres3d1 cresidual=cri;
run;

proc lifetest data=csres3d1 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;
title;

title '3d grade2';
proc lifereg data=oc;
	ods select ModelInfo ParameterEstimates ProbPlot;
	where grade=2;
	model time*event(0) = grade / dist=weibull;
	probplot;
	output out=csres3d2 cresidual=cri;
run;

proc lifetest data=csres3d1 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;
title;

title '3d grade3';
proc lifereg data=oc;
	ods select ModelInfo ParameterEstimates ProbPlot;
	where grade=3;
	model time*event(0) = grade / dist=weibull;
	probplot;
	output out=csres3d3 cresidual=cri;
run;

proc lifetest data=csres3d1 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;
title;

*************************************;
title '3gh control (reference) vs. trmt';
proc sort data=oc;
	by descending trmt;
run;

proc lifereg data=oc order=data;
	ods select ModelInfo ParameterEstimates ProbPlot;
	model time*event(0) = trmt / dist=weibull;
	probplot;
	output out=csres3gh cresidual=cri;
run;

proc lifetest data=csres3gh plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*event(0);
run;
title;

*******************************************************;

* 4 AFT model with Weibull dist
Use Cox-Snell residuals to check if the fitted model is appropriate;

data mel;
	infile '~/BIOS680/data/melanoma.txt';
	input ptid age sex fam rtime rstatus 
	stime sstatus stage $ monilia mumps ppd
	pha sksd trico;
	if upcase(stage) in ('3A', '3B', '3AB') then stages = 3;
	else stages = 4;
run;

proc format;
	value st
	3 = '1-stage3'
	4 = '2-stage4'
;

************************************************;
title '4 control (reference) vs. trmt';
proc lifereg data=mel outest=parm covout order=formatted;
	ods select ModelInfo ParameterEstimates ProbPlot;
	class stages;
	model stime*sstatus(0) = stages / dist=weibull;
	probplot;
	format stages st.;
	output out=csres4 cresidual=cri;
run;

proc lifetest data=csres4 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*sstatus(0);
run;
title;

* 5;
* Create left censor variables in data;
data mel (keep=gender smonth sstatus lower);
	infile '~/BIOS680/data/mela.dat';
	input ptid age gender fam rtime rstatus 
	smonth sstatus stage $ monilia mumps ppd
	pha sksd trico;
	if (gender ^= -9);
	if smonth=. and sstatus=1 then do;
		smonth=30;
		lower=.;
	end;
	else lower=smonth;
run;

  
* Fit AFT Weibull model to examine gender effect on mortality;
* PPlot (percentile plot for Gender1 vs. Gender2) 
checks if AFT is appropriate for 2 groups of data;
title '5 AFT Weibull model left-censoring';
proc lifereg data=mel;
	ods select ParameterEstimates ProbPlot;
	model (lower, smonth)=gender / dist=weibull;
	pplot;
	output out=csres5 cresidual=cri;
run;

* Cox-snell residual plot
(Cumulative hazard H(t)=-logS(t)) vs. residuals) 
rCi have a unit exponential distribution;
* checks if Weibull AFT fit is correct;

proc lifetest data=csres5 plots=(ls) notable;
	ods select NegLogSurvivalPlot;
	time cri*sstatus(0);
run;
title;


