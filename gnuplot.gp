set term post enh color solid
set output "./Output/Timings_vs_residual_GMRES.ps"
set title sprintf("1. NO Preconditioning\t2. JACOBI Preconditioning\t3. GAUSS SEIDEL Preconditioning")
set xlabel 'Time (s)'
set ylabel 'Final Residual'
set grid x y
plot "./Times_FULL.out" u 3:4 w linespoints lc rgb "blue" pt 4 t 'GMRES Full',"./Times_Restarted.out" u 3:4 w linespoints lt 1 lc rgb "red" pt 4 t 'GMRES Restarted'
