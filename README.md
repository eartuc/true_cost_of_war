# True Cost of War
Calculate welfare impact of big shocks from migration flows.

Please check the excel file "Welfare_Calculation_Simple.xls" for a simple demonstration of the methodology, or run "produce_tables.m" in Matlab to reproduce the main tables with many robustness tests. 

The sufficient statistic requires one parameter and one data point: The migration (semi)elasticity parameter $\theta$ and the change in migration flows after the shock, $\Delta \log{m}^{kl}_t$. The paper explains how to convert the raw welfare measure to income equivalent loss and also discusses caveats.

***Formula***\
The upper bound of welfare change in conflict/shock location $k$ (i.e., the lower bound of the welfare impact) can be expressed with the formula below (under the assumption that welfare does not increase in non-conflict/non-shock location $l \neq k$). 
<h2>

</h2>

<h2>
 &emsp;&emsp;&emsp;&emsp;&emsp;  $\Delta{W}^k_t   \leq  \frac{1}{ \theta}  \left(   -  \Delta \log{m}^{kl}_t \right)$   
</h2>

<h3>

</h3>


***Citation***\
Artuc Gomez-Parra and Onder (2026). "True Cost of War: The Conflict in Eastern Ukraine," [Review of Economics and Statistics (forthcoming).](https://doi.org/10.1162/REST.a.1737)

***Working paper version***\
http://documents.worldbank.org/curated/en/099846210202232755


***Abstract***\
This project proposes a new method to estimate the welfare impact of conflicts and remedy common data constraints in conflict-affected environments. The method first estimates how agents regard spatial welfare differentials by voting with their feet, using pre-conflict data. Then, it infers a lower-bound estimate for the conflict-driven welfare shock from partially observed post-conflict migration patterns. Results for the conflict in Eastern Ukraine between 2014 and 2019 show a large lower-bound welfare loss for Donetsk residents equivalent to 8-32 percent of lifetime income depending on agents' time preference and risk aversion parameters.


***Replication package***\
https://doi.org/10.7910/DVN/VSCJWH


***Authors***\
Erhan Artuc (http://www.artuc.org)    
Nicolas Gomez-Parra    
Harun Onder (https://sites.google.com/view/harunonder/)   
