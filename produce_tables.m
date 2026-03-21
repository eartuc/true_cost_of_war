clear

run welfare_main.m
run welfare_weighted.m
run welfare_distweighted.m

disp('Main welfare tables in the paper')
Table2
Table3

save output/Table2.txt Table2 -ascii
save output/Table3.txt Table3 -ascii

