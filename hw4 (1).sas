* 2) Life-table method;

*time midpoint of interval;
*status 1=death, 0=censored, within interval;
* check work;
data life;
	input time status count;
	datalines; 
0.5 1 240
1.5 1 170
1.5 0 10
2.5 1 174
2.5 0 10
3.5 1 138
4.5 1 98
4.5 0 20
5.5 1 60
6.5 1 52
7.5 1 44
8.5 1 32
9 1 28
9 0 24
;

proc lifetest method=life 
	intervals=(0 to 9 by 1)
	plots=(hazard, survival, pdf);
	time time*status(0);
	freq count;
run;
*/
* 3d-i;

proc import out=oc
    datafile="~/BIOS680/data/oc.txt"
    dbms=dlm
    replace;
    getnames=no;
    delimiter=' ';
    guessingrows=max;
run;

data surv;
	set oc;
	rename var5=group;
	rename var6=time;
	rename var7=delta;
	rename var3=stage;
	label	
		time = 'Time from second look laparotomy to recurrence (days)';
		delta = 'Event indicator: 0=censored, 1=treat';
		group = 'Treatment: 0=not given, 1 = observed';
run;

* 3d Nelson estimator for cumulative hazard;
* output estimates;
ods graphics on;
ods output ProductLimitEstimates = est;
proc lifetest data=surv outs=outsurv1 noprint
	nelson alpha=0.05;
	strata group;
	time time * delta(0);
run;

* 3d Plot Nelson cumulative hazard;
title '3d. Nelson cumulative hazard estimate for 2 groups from OC data';
proc sgplot data = est;
	series x = time y = CumHaz / group=group;
run;

* 3e Plot kernel hazard function estimates 
using bandwidth for each treatment group;

* bandwidth 100;
title '3e. bandwidth 100';
ods graphics on;
ods select HazardPlot;
proc lifetest data=surv 
	plots=hazard(BW=100);
	time time * delta(0);
	strata group;
run;
title;

* 3e;
* Kernel-hazard function bandwidth of 100, 40, 130;
ods select HazardPlot;
title '3e. bandwidth 40';
proc lifetest data=surv 
	plots=hazard(BW=40);
	time time * delta(0);
	strata group;
run;
title;

ods select HazardPlot;
title '3e. bandwidth 130';
proc lifetest data=surv 
	plots=hazard(BW=130);
	time time * delta(0);
	strata group;
run;
title;

* 3h;
* Logrank test;
title '3h. logrank test';
ods select HomTests SurvivalPlot;
proc lifetest data=surv plots=s(test);
	strata group / test=logrank;
	time time * delta(0);
run;
title;

* 3i;
* Gehan/Breslow generalized Wilcoxon test;
* shows all tests by default;
title '3i. logrank test';
ods select HomTests SurvivalPlot;
proc lifetest data=surv plots=s(test);
	strata group;
	time time * delta(0);
run;
title;

*4a-c; 
title '4a. Kaplan-Meier survival curves';
ods select SurvivalPlot;
proc lifetest data=surv
	plots=survival(cl cb=all);
	strata group ;
	time time * delta(0);
	survival out=loglog conftype=loglog;
run;
title;

*4d;
* Cumulative hazard for 3 stage groups;
*ods output ProductLimitEstimates=est2; *displays pointwise CI;
proc lifetest data=surv outs=outsurv2 plots=survival(cl cb=all) 
	nelson alpha=0.05;
	strata stage;
	time time * delta(0);
run;

* Plot Nelson cumulative hazard;
title '4d. Nelson cumulative hazard estimate for 3 stages from OC data';
proc sgplot data = est2;
	series x = time y = CumHaz / group=stage;
run;
title;

ods graphics off;