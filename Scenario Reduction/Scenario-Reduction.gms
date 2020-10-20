****Multi-Stage Farmer Problem-Using Scenario Reduction-
set crop/wheat, corn, sugarbeet/;
set scenario/1*20/;
alias(scenario,scenario1,scenario2);

$CALL GDXXRW.EXE probsr.xls par=prob rng=a1:b20 rdim=1
parameter prob(scenario);
$gdxin probsr.gdx
$load prob
$gdxin

$CALL GDXXRW.EXE yield.xls par=yield rng=a1:u4
parameter yield(crop,scenario);
$gdxin yield.gdx
$load yield
$gdxin

parameter plantcost(crop)/wheat 150, corn 230, sugarbeet 260/;
parameter saleprice(crop)/wheat 170, corn 150, sugarbeet 36/;
parameter xsellprice(crop)/sugarbeet 10/;
parameter purchaseprice(crop)/wheat 238, corn 210/;
parameter minreq(crop)/wheat 200, corn 240/;
scalar land/500/;
scalar maxbeet/6000/;
****
*
****

variables
surface1(crop), sales1(crop,scenario1),xsales1(crop,scenario1),
purchase1(crop,scenario1),surface2(crop,scenario1), sales2(crop,scenario1,scenario2),xsales2(crop,scenario1,scenario2),
purchase2(crop,scenario1,scenario2),totcost;
positive variables
surface1(crop), sales1(crop,scenario1),xsales1(crop,scenario1),purchase1(crop,scenario1),
surface2(crop,scenario1), sales2(crop,scenario1,scenario2),xsales2(crop,scenario1,scenario2),
purchase2(crop,scenario1,scenario2);
free variable totcost;
sales1.up('sugarbeet',scenario1)=maxbeet;
sales2.up('sugarbeet',scenario1,scenario2)=maxbeet;

equations
totobj, maxland1, minlvl1(crop,scenario1), maxland2(scenario1), minlvl2(crop,scenario1,scenario2),trenial(scenario1);

totobj..totcost=e=sum(crop,plantcost(crop)*surface1(crop))
  -sum((crop,scenario1),prob(scenario1)*saleprice(crop)*sales1(crop,scenario1))-sum((crop,scenario1) $ xsellprice(crop),prob(scenario1)*xsellprice(crop)*xsales1(crop,scenario1))
  +sum((crop,scenario1) $ purchaseprice(crop),prob(scenario1)*purchaseprice(crop)*purchase1(crop,scenario1))
  +sum((crop,scenario1),prob(scenario1)*surface2(crop,scenario1))
  -sum((crop,scenario1,scenario2),prob(scenario1)*prob(scenario2)*saleprice(crop)*sales2(crop,scenario1,scenario2))-sum((crop,scenario1,scenario2) $ xsellprice(crop),prob(scenario1)*prob(scenario2)*xsellprice(crop)*xsales2(crop,scenario1,scenario2))
  +sum((crop,scenario1,scenario2) $ purchaseprice(crop),prob(scenario1)*prob(scenario2)*purchaseprice(crop)*purchase2(crop,scenario1,scenario2))
;

maxland1..sum(crop,surface1(crop))=l=land;
minlvl1(crop,scenario1)..yield(crop,scenario1)*surface1(crop)+purchase1(crop,scenario1)$purchaseprice(crop)=g=
              sales1(crop,scenario1)+xsales1(crop,scenario1)$xsellprice(crop)+minreq(crop);
trenial(scenario1)..surface2('sugarbeet',scenario1)=l=500-surface1('sugarbeet');
maxland2(scenario1)..sum(crop,surface2(crop,scenario1))=l=land;
minlvl2(crop,scenario1,scenario2)..yield(crop,scenario1)*surface2(crop,scenario1)+purchase2(crop,scenario1,scenario2)$purchaseprice(crop)=g=
              sales2(crop,scenario1,scenario2)+xsales2(crop,scenario1,scenario2)$xsellprice(crop)+minreq(crop);



model mainproblem/all/;
option optcr=0;
option lp=cplex;
solve mainproblem using lp minimizing totcost;
display totcost.l,surface1.l,surface2.l,sales1.l,xsales1.l,purchase1.l,sales2.l,xsales2.l,purchase2.l;