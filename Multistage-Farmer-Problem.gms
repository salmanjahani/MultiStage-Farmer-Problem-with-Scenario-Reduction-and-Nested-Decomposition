****Multi-Stage Farmer Problem****Final***
set iter loopiteration/1*20/
    cnew1(iter) Dynamic set
    cnew21(iter) Dynamic set
    cnew22(iter) Dynamic set
    cnew23(iter) Dynamic set;
    cnew1(iter)=no;
    cnew21(iter)=no;
    cnew22(iter)=no;
    cnew23(iter)=no;
set crop/wheat, corn, sugarbeet/;
set scenario/above, normal, below/;
alias(scenario,scenario2,scenario3);
set scenario1(scenario) /normal/;
set period stage/1,2,3/;
sets
eq2 /1*6/
eq3 /1*4/;
parameter meanyield(crop)/wheat 2.5, corn 3, sugarbeet 20/;
parameter fraction(scenario)/above 0.2, normal 0, below -0.2/;
parameter yield(crop,scenario);
yield(crop,scenario)=(1+fraction(scenario))*meanyield(crop);
parameter prob(scenario);
prob(scenario)=1/3;
parameter plantcost(crop)/wheat 150, corn 230, sugarbeet 260/;
parameter saleprice(crop)/wheat 170, corn 150, sugarbeet 36/;
parameter xsellprice(crop)/sugarbeet 10/;
parameter purchaseprice(crop)/wheat 238, corn 210/;
parameter minreq(crop)/wheat 200, corn 240/;
scalar land/500/;
scalar maxbeet/6000/;
parameter yieldsc(crop);
parameter surfacesc(crop);
********************
**Cuts' Parameters**
********************
parameters e1(crop,iter),e2b(crop,iter),e2n(crop,iter),e2a(crop,iter);
parameters e11(iter),e22b(iter),e22n(iter),e22a(iter);

variables
surface(crop), sales(crop),xsales(crop),
purchase(crop),totcost,teta;
positive variables
surface, sales, xsales, purchase,teta;
free variable totcost;

equations
totobjm1,totobj1,
totobjm,totobj,totobjl,
maxland, minlvl(crop),beet,trenial,
cut1(iter),cut2b(iter),cut2n(iter),cut2a(iter);

totobjm1..totcost=e=sum(crop,plantcost(crop)*surface(crop));
totobj1..totcost=e=sum(crop,plantcost(crop)*surface(crop))
  -sum(crop,saleprice(crop)*sales(crop))-sum(crop $ xsellprice(crop),xsellprice(crop)*xsales(crop))
  +sum(crop $ purchaseprice(crop),purchaseprice(crop)*purchase(crop));

totobjm..totcost=e=sum(crop,plantcost(crop)*surface(crop))+teta;
totobj..totcost=e=sum(crop,plantcost(crop)*surface(crop))
  -sum(crop,saleprice(crop)*sales(crop))-sum(crop $ xsellprice(crop),xsellprice(crop)*xsales(crop))
  +sum(crop $ purchaseprice(crop),purchaseprice(crop)*purchase(crop))+teta;

totobjl..totcost=e=
  -sum(crop,saleprice(crop)*sales(crop))-sum(crop $ xsellprice(crop),xsellprice(crop)*xsales(crop))
  +sum(crop $ purchaseprice(crop),purchaseprice(crop)*purchase(crop));

maxland..sum(crop,surface(crop))=l=land;
minlvl(crop)..yieldsc(crop)*surfacesc(crop)+purchase(crop)$purchaseprice(crop)=g=
              sales(crop)+xsales(crop)$xsellprice(crop)+minreq(crop);
beet..sales('sugarbeet')=l=maxbeet;
trenial..surface('sugarbeet')+surfacesc('sugarbeet')=l=500;

cut2b(cnew21)..sum(crop,e2b(crop,cnew21)*surface(crop))+teta=g=e22b(cnew21);
cut2n(cnew22)..sum(crop,e2n(crop,cnew22)*surface(crop))+teta=g=e22n(cnew22);
cut2a(cnew23)..sum(crop,e2a(crop,cnew23)*surface(crop))+teta=g=e22a(cnew23);

cut1(cnew1)..sum(crop,e1(crop,cnew1)*surface(crop))+teta=g=e11(cnew1);

model masterproblem1/totobjm1, maxland/;
model subproblem1/totobj1,maxland,minlvl,beet,trenial/;

model masterproblem/totobjm, maxland,cut1/;
model subproblemb/totobj,maxland,minlvl,beet,trenial,cut2b/;
model subproblemn/totobj,maxland,minlvl,beet,trenial,cut2n/;
model subproblema/totobj,maxland,minlvl,beet,trenial,cut2a/;

model finalpro/totobjl,minlvl,beet/;
option optcr=0;
option lp=baron;
****************************
****Nested Decomposition****
****************************
scalar converged /0/;
parameters pi3(iter,scenario2,scenario3,eq3),pi2(iter,scenario2,eq2);
parameters sigmab(iter),sigman(iter),sigmaa(iter);
parameters surfacep1(iter,crop),surfacep2(iter,crop,scenario2),salesp(iter,period,crop,scenario),
xsalesp(iter,period,crop,scenario),purchasep(iter,period,crop,scenario),tetap1,tetap2(scenario2);
parameters tetabar11, tetabar22b,tetabar22n,tetabar22a;
****Initialization****
solve masterproblem1 using lp minimizing totcost;
surfacep1('1',crop)=surface.l(crop);
tetap1=0;
 loop(scenario2,
  yieldsc(crop)=yield(crop,scenario2);
  surfacesc(crop)=surfacep1('1',crop);
  solve subproblem1 using lp minimizing totcost;
  surfacep2('1',crop,scenario2)=surface.l(crop);
  pi2('1',scenario2,'1')=maxland.m;
  pi2('1',scenario2,'2')=minlvl.m('wheat');
  pi2('1',scenario2,'3')=minlvl.m('corn');
  pi2('1',scenario2,'4')=minlvl.m('sugarbeet');
  pi2('1',scenario2,'5')=beet.m;
  pi2('1',scenario2,'6')=trenial.m;
  tetap2(scenario2)=0;
 );
loop(scenario2,
 loop(scenario3,
 surfacesc(crop)=surfacep2('1',crop,scenario2);
 yieldsc(crop)=yield(crop,scenario3);
 solve finalpro using lp minimizing totcost;
 pi3('1',scenario2,scenario3,'1')=minlvl.m('wheat');
 pi3('1',scenario2,scenario3,'2')=minlvl.m('corn');
 pi3('1',scenario2,scenario3,'3')=minlvl.m('sugarbeet');
 pi3('1',scenario2,scenario3,'4')=beet.m;
 );
);
****OK****
display surfacep1,surfacep2,pi2,pi3;
*********************************
****Begining of Backward Pass****
*********Cut Computation*********
*********************************
  e2b('wheat','1')=sum(scenario3,prob(scenario3)*yield('wheat',scenario3)*pi3('1','below',scenario3,'1'));
  e2n('wheat','1')=sum(scenario3,prob(scenario3)*yield('wheat',scenario3)*pi3('1','normal',scenario3,'1'));
  e2a('wheat','1')=sum(scenario3,prob(scenario3)*yield('wheat',scenario3)*pi3('1','above',scenario3,'1'));
  e2b('corn','1')=sum(scenario3,prob(scenario3)*yield('corn',scenario3)*pi3('1','below',scenario3,'2'));
  e2n('corn','1')=sum(scenario3,prob(scenario3)*yield('corn',scenario3)*pi3('1','normal',scenario3,'2'));
  e2a('corn','1')=sum(scenario3,prob(scenario3)*yield('corn',scenario3)*pi3('1','above',scenario3,'2'));
  e2b('sugarbeet','1')=sum(scenario3,prob(scenario3)*yield('sugarbeet',scenario3)*pi3('1','below',scenario3,'3'));
  e2n('sugarbeet','1')=sum(scenario3,prob(scenario3)*yield('sugarbeet',scenario3)*pi3('1','normal',scenario3,'3'));
  e2a('sugarbeet','1')=sum(scenario3,prob(scenario3)*yield('sugarbeet',scenario3)*pi3('1','above',scenario3,'3'));
  e22b('1')=sum(scenario3,prob(scenario3)*(200*pi3('1','below',scenario3,'1')+240*pi3('1','below',scenario3,'2')+6000*pi3('1','below',scenario3,'4')));
  e22n('1')=sum(scenario3,prob(scenario3)*(200*pi3('1','normal',scenario3,'1')+240*pi3('1','normal',scenario3,'2')+6000*pi3('1','normal',scenario3,'4')));
  e22a('1')=sum(scenario3,prob(scenario3)*(200*pi3('1','above',scenario3,'1')+240*pi3('1','above',scenario3,'2')+6000*pi3('1','above',scenario3,'4')));
  tetabar22b=e22b('1')-sum(crop,e2b(crop,'1')*surfacep2('1',crop,'below'));
  tetabar22n=e22n('1')-sum(crop,e2n(crop,'1')*surfacep2('1',crop,'normal'));
  tetabar22a=e22a('1')-sum(crop,e2a(crop,'1')*surfacep2('1',crop,'above'));
****
  if(tetabar22b>tetap2('below'),
     cnew21('1')=yes;
      yieldsc(crop)=yield(crop,'below');
      surfacesc(crop)=surfacep1('1',crop);
      solve subproblemb using lp minimizing totcost;
      surfacep2('1',crop,'below')=surface.l(crop);
      pi2('1','below','1')=maxland.m;
      pi2('1','below','2')=minlvl.m('wheat');
      pi2('1','below','3')=minlvl.m('corn');
      pi2('1','below','4')=minlvl.m('sugarbeet');
      pi2('1','below','5')=beet.m;
      pi2('1','below','6')=trenial.m;
      sigmab('1')=cut2b.m('1');
   );

  if(tetabar22n>tetap2('normal'),
   cnew22('1')=yes;
      yieldsc(crop)=yield(crop,'normal');
      surfacesc(crop)=surfacep1('1',crop);
      solve subproblemn using lp minimizing totcost;
      surfacep2('1',crop,'normal')=surface.l(crop);
      pi2('1','normal','1')=maxland.m;
      pi2('1','normal','2')=minlvl.m('wheat');
      pi2('1','normal','3')=minlvl.m('corn');
      pi2('1','normal','4')=minlvl.m('sugarbeet');
      pi2('1','normal','5')=beet.m;
      pi2('1','normal','6')=trenial.m;
      sigman('1')=cut2n.m('1');
  );
  if(tetabar22a>tetap2('above'),
    cnew23('1')=yes;
      yieldsc(crop)=yield(crop,'above');
      surfacesc(crop)=surfacep1('1',crop);
      solve subproblema using lp minimizing totcost;
      surfacep2('1',crop,'above')=surface.l(crop);
      pi2('1','above','1')=maxland.m;
      pi2('1','above','2')=minlvl.m('wheat');
      pi2('1','above','3')=minlvl.m('corn');
      pi2('1','above','4')=minlvl.m('sugarbeet');
      pi2('1','above','5')=beet.m;
      pi2('1','above','6')=trenial.m;
      sigmaa('1')=cut2a.m('1');
  );

****Add Cut****
****Recalculation****
****End of Stage II****
    e1('wheat','1')=sum(scenario2,prob(scenario2)*yield('wheat',scenario2)*pi2('1',scenario2,'2'));
    e1('corn','1')=sum(scenario2,prob(scenario2)*yield('corn',scenario2)*pi2('1',scenario2,'3'));
    e1('sugarbeet','1')=sum(scenario2,prob(scenario2)*(yield('sugarbeet',scenario2)*pi2('1',scenario2,'4')+pi2('1',scenario2,'6')));
    e11('1')=sum(scenario2,prob(scenario2)*(500*pi2('1',scenario2,'1')+200*pi2('1',scenario2,'2')
                  +240*pi2('1',scenario2,'3')+6000*pi2('1',scenario2,'5')+500*pi2('1',scenario2,'6')))
                  +prob('below')*(sum(cnew21,sigmab(cnew21)*e22b(cnew21)))
                  +prob('normal')*(sum(cnew22,sigman(cnew22)*e22n(cnew22)))
                  +prob('above')*(sum(cnew23,sigmaa(cnew23)*e22a(cnew23)))  ;
    tetabar11=e11('1')-sum(crop,e1(crop,'1')*surfacep1('1',crop));

****End of Stage I****
    if(tetabar11>tetap1,
      cnew1('1')=yes;
    else
      converged=1;
      display surfacep1,surfacep2,converged,tetap1,tetap2;
      execute_unload "result.gdx" surfacep1,converged;
      execute 'gdxxrw.exe result.gdx par=surfacep rng=Surface!'
      execute 'gdxxrw.exe result.gdx par=converged rng=Converged!'
      Abort$(converged=1) 'Algorithm has converged';
     );
*****************************
******End of Backward Pass***
*****************************
****Main Loop Computation****
*****************************
loop(iter $ (ord(iter)> 1 and converged ne 1),
********************************
****Begining of Forward Pass****
********************************
****Stage I****
 solve masterproblem using lp minimizing totcost;
 surfacep1(iter,crop)=surface.l(crop);
 tetap1=teta.l;
****Stage II****
  yieldsc(crop)=yield(crop,'below');
  surfacesc(crop)=surfacep1(iter,crop);
  solve subproblemb using lp minimizing totcost;
  surfacep2(iter,crop,'below')=surface.l(crop);
  pi2(iter,'below','1')=maxland.m;
  pi2(iter,'below','2')=minlvl.m('wheat');
  pi2(iter,'below','3')=minlvl.m('corn');
  pi2(iter,'below','4')=minlvl.m('sugarbeet');
  pi2(iter,'below','5')=beet.m;
  pi2(iter,'below','6')=trenial.m;
  tetap2('below')=teta.l;

  yieldsc(crop)=yield(crop,'normal');
  surfacesc(crop)=surfacep1(iter,crop);
  solve subproblemn using lp minimizing totcost;
  surfacep2(iter,crop,'normal')=surface.l(crop);
  pi2(iter,'normal','1')=maxland.m;
  pi2(iter,'normal','2')=minlvl.m('wheat');
  pi2(iter,'normal','3')=minlvl.m('corn');
  pi2(iter,'normal','4')=minlvl.m('sugarbeet');
  pi2(iter,'normal','5')=beet.m;
  pi2(iter,'normal','6')=trenial.m;
  tetap2('normal')=teta.l;

  yieldsc(crop)=yield(crop,'above');
  surfacesc(crop)=surfacep1(iter,crop);
  solve subproblema using lp minimizing totcost;
  surfacep2(iter,crop,'above')=surface.l(crop);
  pi2(iter,'above','1')=maxland.m;
  pi2(iter,'above','2')=minlvl.m('wheat');
  pi2(iter,'above','3')=minlvl.m('corn');
  pi2(iter,'above','4')=minlvl.m('sugarbeet');
  pi2(iter,'above','5')=beet.m;
  pi2(iter,'above','6')=trenial.m;
  tetap2('above')=teta.l;

****Stage III****
loop(scenario2,
 loop(scenario3,
 yieldsc(crop)=yield(crop,scenario3);
 surfacesc(crop)=surfacep2(iter,crop,scenario2);
 solve finalpro using lp minimizing totcost;
 pi3(iter,scenario2,scenario3,'1')=minlvl.m('wheat');
 pi3(iter,scenario2,scenario3,'2')=minlvl.m('corn');
 pi3(iter,scenario2,scenario3,'3')=minlvl.m('sugarbeet');
 pi3(iter,scenario2,scenario3,'4')=beet.m;
 );
);
***************************
****End of Forward Pass****
***************************
*********************************
****Begining of Backward Pass****
*********Cut Computation*********
*********************************
  e2b('wheat',iter)=sum(scenario3,prob(scenario3)*yield('wheat',scenario3)*pi3(iter,'below',scenario3,'1'));
  e2n('wheat',iter)=sum(scenario3,prob(scenario3)*yield('wheat',scenario3)*pi3(iter,'normal',scenario3,'1'));
  e2a('wheat',iter)=sum(scenario3,prob(scenario3)*yield('wheat',scenario3)*pi3(iter,'above',scenario3,'1'));
  e2b('corn',iter)=sum(scenario3,prob(scenario3)*yield('corn',scenario3)*pi3(iter,'below',scenario3,'2'));
  e2n('corn',iter)=sum(scenario3,prob(scenario3)*yield('corn',scenario3)*pi3(iter,'normal',scenario3,'2'));
  e2a('corn',iter)=sum(scenario3,prob(scenario3)*yield('corn',scenario3)*pi3(iter,'above',scenario3,'2'));
  e2b('sugarbeet',iter)=sum(scenario3,prob(scenario3)*yield('sugarbeet',scenario3)*pi3(iter,'below',scenario3,'3'));
  e2n('sugarbeet',iter)=sum(scenario3,prob(scenario3)*yield('sugarbeet',scenario3)*pi3(iter,'normal',scenario3,'3'));
  e2a('sugarbeet',iter)=sum(scenario3,prob(scenario3)*yield('sugarbeet',scenario3)*pi3(iter,'above',scenario3,'3'));
  e22b(iter)=sum(scenario3,prob(scenario3)*(200*pi3(iter,'below',scenario3,'1')+240*pi3(iter,'below',scenario3,'2')+6000*pi3(iter,'below',scenario3,'4')));
  e22n(iter)=sum(scenario3,prob(scenario3)*(200*pi3(iter,'normal',scenario3,'1')+240*pi3(iter,'normal',scenario3,'2')+6000*pi3(iter,'normal',scenario3,'4')));
  e22a(iter)=sum(scenario3,prob(scenario3)*(200*pi3(iter,'above',scenario3,'1')+240*pi3(iter,'above',scenario3,'2')+6000*pi3(iter,'above',scenario3,'4')));
  tetabar22b=e22b(iter)-sum(crop,e2b(crop,iter)*surfacep2(iter,crop,'below'));
  tetabar22n=e22n(iter)-sum(crop,e2n(crop,iter)*surfacep2(iter,crop,'normal'));
  tetabar22a=e22a(iter)-sum(crop,e2a(crop,iter)*surfacep2(iter,crop,'above'));
****
  if(tetabar22b>tetap2('below'),
     cnew21(iter)=yes;
      yieldsc(crop)=yield(crop,'below');
      surfacesc(crop)=surfacep1(iter,crop);
      solve subproblemb using lp minimizing totcost;
      surfacep2(iter,crop,'below')=surface.l(crop);
      pi2(iter,'below','1')=maxland.m;
      pi2(iter,'below','2')=minlvl.m('wheat');
      pi2(iter,'below','3')=minlvl.m('corn');
      pi2(iter,'below','4')=minlvl.m('sugarbeet');
      pi2(iter,'below','5')=beet.m;
      pi2(iter,'below','6')=trenial.m;
      sigmab(iter)=cut2b.m(iter);
  );

  if(tetabar22n>tetap2('normal'),
   cnew22(iter)=yes;
      yieldsc(crop)=yield(crop,'normal');
      surfacesc(crop)=surfacep1(iter,crop);
      solve subproblemn using lp minimizing totcost;
      surfacep2(iter,crop,'normal')=surface.l(crop);
      pi2(iter,'normal','1')=maxland.m;
      pi2(iter,'normal','2')=minlvl.m('wheat');
      pi2(iter,'normal','3')=minlvl.m('corn');
      pi2(iter,'normal','4')=minlvl.m('sugarbeet');
      pi2(iter,'normal','5')=beet.m;
      pi2(iter,'normal','6')=trenial.m;
      sigman(iter)=cut2n.m(iter);
  );
  if(tetabar22a>tetap2('above'),
    cnew23(iter)=yes;
      yieldsc(crop)=yield(crop,'above');
      surfacesc(crop)=surfacep1(iter-1,crop);
      solve subproblema using lp minimizing totcost;
      surfacep2(iter,crop,'above')=surface.l(crop);
      pi2(iter,'above','1')=maxland.m;
      pi2(iter,'above','2')=minlvl.m('wheat');
      pi2(iter,'above','3')=minlvl.m('corn');
      pi2(iter,'above','4')=minlvl.m('sugarbeet');
      pi2(iter,'above','5')=beet.m;
      pi2(iter,'above','6')=trenial.m;
      sigmaa(iter)=cut2a.m(iter);
  );

****Add Cut****
****Recalculation****
****End of Stage II****
    e1('wheat',iter)=sum(scenario2,prob(scenario2)*yield('wheat',scenario2)*pi2(iter,scenario2,'2'));
    e1('corn',iter)=sum(scenario2,prob(scenario2)*yield('corn',scenario2)*pi2(iter,scenario2,'3'));
    e1('sugarbeet',iter)=sum(scenario2,prob(scenario2)*(yield('sugarbeet',scenario2)*pi2(iter,scenario2,'4')+pi2(iter,scenario2,'6')));
    e11(iter)=sum(scenario2,prob(scenario2)*(500*pi2(iter,scenario2,'1')+200*pi2(iter,scenario2,'2')
                  +240*pi2(iter,scenario2,'3')+6000*pi2(iter,scenario2,'5')+500*pi2(iter,scenario2,'6')))
                  +prob('below')*(sum(cnew21,sigmab(cnew21)*e22b(cnew21)))
                  +prob('normal')*(sum(cnew22,sigman(cnew22)*e22n(cnew22)))
                  +prob('above')*(sum(cnew23,sigmaa(cnew23)*e22a(cnew23)))  ;
    tetabar11=e11(iter)-sum(crop,e1(crop,iter)*surfacep1(iter,crop));

****End of Stage I****
    if(tetabar11>tetap1,
      cnew1(iter)=yes;
    else
      converged=1;
      display surfacep1,surfacep2,converged,tetap1,tetap2;
      execute_unload "result.gdx" surfacep1,converged;
      execute 'gdxxrw.exe result.gdx par=surfacep rng=Surface!'
      execute 'gdxxrw.exe result.gdx par=converged rng=Converged!'
      Abort$(converged=1) 'Algorithm has converged';
     );
*****************************
******End of Backward Pass***
*****************************
);
*********************************************
*************End of Main Loop****************
****End of Nested Decomposition Algorithm****
*********************************************