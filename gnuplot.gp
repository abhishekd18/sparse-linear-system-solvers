set term post enh color solid
set output "./Output/Timings_GMRES.ps"
set title sprintf("1. NO Preconditioning\t2. JACOBI Preconditioning\t3. GAUSS SEIDEL Preconditioning")
set xlabel 'Preconditioning'
set ylabel 'Time (s)'
set grid x y
set xtics 1,1,3
set key right
plot "./Times_FULL.out" u 1:3 w linespoints lc rgb "blue" pt 4 t 'GMRES Full',"./Times_Restarted.out" u 1:3 w linespoints lt 1 lc rgb "red" pt 4 t 'GMRES Restarted'
