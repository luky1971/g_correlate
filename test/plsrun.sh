g_correlate -f pdz50ns_nowat.xtc -s pdz50ns.tpr -a 'N H'
echo 'Diffing unit vectors'
diff vecs_fbyp.txt NHvecs_fbyp.txt
diff vecs_pbyf.txt NHvecs_pbyf.txt
