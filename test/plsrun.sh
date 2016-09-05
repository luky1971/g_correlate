g_correlate -f pdz50ns_nowat.xtc -s pdz50ns.tpr -a 'N H' -nt 1000
#echo 'Diffing unit vectors'
#diff vecs_fbyp.txt NHvecs_fbyp.txt
#diff vecs_pbyf.txt NHvecs_pbyf.txt
echo 'Diffing autocorrelations'
diff N-H_corr.dat N-Hexp_corr_nt1000.dat 
echo 'Diffing S2s'
diff s2.dat N-Hexp_s2_nt1000.dat
