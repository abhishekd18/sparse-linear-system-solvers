rm ./Output/*Error*
rm *.out
./bin/main orsirr_1.mtx GMRES_FULL NULL>run_gmres_full_no.out
echo "set term post enh color solid
set output \"./Output/GMRES_FULL_NO.ps\"
set title 'GMRES FULL No Preconditioning'
set xlabel 'Iterations'
set ylabel 'Relative Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_FULL_NO.ps ./Output/GMRES_FULL_NO.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_FULL_NO.out

./bin/main orsirr_1.mtx GMRES_FULL JACOBI>run_gmres_full_jacobi.out
echo "set term post enh color solid
set output \"./Output/GMRES_FULL_JACOBI.ps\"
set title 'GMRES FULL JACOBI'
set xlabel 'Iterations'
set ylabel 'Relative Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_FULL_JACOBI.ps ./Output/GMRES_FULL_JACOBI.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_FULL_JACOBI.out

./bin/main orsirr_1.mtx GMRES_FULL GAUSS_SEIDEL>run_gmres_full_gauss_seidel.out
echo "set term post enh color solid
set output \"./Output/GMRES_FULL_GAUSS_SEIDEL.ps\"
set title 'GMRES FULL GAUSS SEIDEL'
set xlabel 'Iterations'
set ylabel 'Relative Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lt 1 lc rgb \"blue\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_FULL_GAUSS_SEIDEL.ps ./Output/GMRES_FULL_GAUSS_SEIDEL.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_FULL_GAUSS_SEIDEL.out

./bin/main orsirr_1.mtx GMRES_RESTARTED NULL>run_gmres_restarted_no.out
echo "set term post enh color solid
set output \"./Output/GMRES_RESTARTED_NO.ps\"
set title 'GMRES RESTARTED No Preconditioning'
set xlabel 'Iterations'
set ylabel 'Absolute Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lw 0.5 lt 1 lc rgb \"red\" title \"Relative Residual-Internal Iterations\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lw 0.5 lt 1 lc rgb \"blue\" title \"Absolute Residual-Internal Iterations\",\
\"./Output/Error_GMRES_Restarted.out\" u 1:3 w l axes x1y1 lw 2 lt 1 lc rgb \"#8B0000\" title \"Relative Residual\",\
\"./Output/Error_GMRES_Restarted.out\" u 1:4 w l axes x1y1 lw 2 lt 1 lc rgb \"#00008B\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_RESTARTED_NO.ps ./Output/GMRES_RESTARTED_NO.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_RESTARTED_NO_internal.out
mv ./Output/Error_GMRES_Restarted.out ./Output/Error_GMRES_RESTARTED_NO.out

./bin/main orsirr_1.mtx GMRES_RESTARTED JACOBI>run_gmres_restarted_jacobi.out
echo "set term post enh color solid
set output \"./Output/GMRES_RESTARTED_JACOBI.ps\"
set title 'GMRES RESTARTED JACOBI'
set xlabel 'Iterations'
set ylabel 'Absolute Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lw 0.5 lt 1 lc rgb \"red\" title \"Relative Residual-Internal Iterations\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lw 0.5 lt 1 lc rgb \"blue\" title \"Absolute Residual-Internal Iterations\",\
\"./Output/Error_GMRES_Restarted.out\" u 1:3 w l axes x1y1 lw 2 lt 1 lc rgb \"#8B0000\" title \"Relative Residual\",\
\"./Output/Error_GMRES_Restarted.out\" u 1:4 w l axes x1y1 lw 2 lt 1 lc rgb \"#00008B\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_RESTARTED_JACOBI.ps ./Output/GMRES_RESTARTED_JACOBI.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_RESTARTED_JACOBI_internal.out
mv ./Output/Error_GMRES_Restarted.out ./Output/Error_GMRES_RESTARTED_JACOBI.out

./bin/main orsirr_1.mtx GMRES_RESTARTED GAUSS_SEIDEL>run_gmres_restarted_gauss_seidel.out
echo "set term post enh color solid
set output \"./Output/GMRES_RESTARTED_GAUSS_SEIDEL.ps\"
set title 'GMRES RESTARTED GAUSS SEIDEL'
set xlabel 'Iterations'
set ylabel 'Absolute Residual'
set y2label 'Error ||x* - xm||'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_GMRES.out\" u 1:2 w l axes x1y2 lt 1 lc rgb \"green\" title \"Error ||x* - xm||\",\
\"./Output/Error_GMRES.out\" u 1:3 w l axes x1y1 lw 0.5 lt 1 lc rgb \"red\" title \"Relative Residual-Internal Iterations\",\
\"./Output/Error_GMRES.out\" u 1:4 w l axes x1y1 lw 0.5 lt 1 lc rgb \"blue\" title \"Absolute Residual-Internal Iterations\",\
\"./Output/Error_GMRES_Restarted.out\" u 1:3 w l axes x1y1 lw 2 lt 1 lc rgb \"#8B0000\" title \"Relative Residual\",\
\"./Output/Error_GMRES_Restarted.out\" u 1:4 w l axes x1y1 lw 2 lt 1 lc rgb \"#00008B\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/GMRES_RESTARTED_GAUSS_SEIDEL.ps ./Output/GMRES_RESTARTED_GAUSS_SEIDEL.pdf
mv ./Output/Error_GMRES.out ./Output/Error_GMRES_RESTARTED_GAUSS_SEIDEL_internal.out
mv ./Output/Error_GMRES_Restarted.out ./Output/Error_GMRES_RESTARTED_GAUSS_SEIDEL.out

./bin/main s3rmt3m3.mtx CG>run_cg.out
echo "set term post enh color solid
set output \"./Output/CG.ps\"
set title 'Conjugate Gradient'
set xlabel 'Iterations'
set ylabel 'Relative Residual'
set y2label 'Error'
set y2tics
set log y
set log y2
set grid x y
plot \"./Output/Error_CG.out\" u 1:2 w l axes x1y2 lw 2 lt 1 lc rgb \"green\" title \"Error in A norm\",\
\"./Output/Error_CG.out\" u 1:3 w l axes x1y2 lw 2 lt 1 lc rgb \"blue\" title \"Error in 2 norm\",\
\"./Output/Error_CG.out\" u 1:4 w l axes x1y1 lw 0.1 lt 1 lc rgb \"red\" title \"Relative Residual\",\
\"./Output/Error_CG.out\" u 1:5 w l axes x1y1 lw 0.1 lt 1 lc rgb \"magenta\" title \"Absolute Residual\"">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/CG.ps ./Output/CG.pdf

grep -e "Time" "./run_gmres_full_no.out"|cut -d"=" -f4|cut -d"s" -f1>>timings.out
grep -e "Time" "./run_gmres_full_jacobi.out"|cut -d"=" -f4|cut -d"s" -f1>>timings.out
grep -e "Time" "./run_gmres_full_gauss_seidel.out"|cut -d"=" -f4|cut -d"s" -f1>>timings.out
paste preconditioners timings.out>Times_FULL.out
grep -e "Time" "./run_gmres_restarted_no.out"|cut -d"=" -f4|cut -d"s" -f1>>timings_res.out
grep -e "Time" "./run_gmres_restarted_jacobi.out"|cut -d"=" -f4|cut -d"s" -f1>>timings_res.out
grep -e "Time" "./run_gmres_restarted_gauss_seidel.out"|cut -d"=" -f4|cut -d"s" -f1>>timings_res.out
paste preconditioners timings_res.out>Times_Restarted.out
echo "set term post enh color solid
set output \"./Output/Timings_GMRES.ps\"
set title sprintf(\"1. NO Preconditioning\t2. JACOBI Preconditioning\t3. GAUSS SEIDEL Preconditioning\")
set xlabel 'Preconditioning'
set ylabel 'Time (s)'
set grid x y
set xtics 1,1,3
set key right
plot \"./Times_FULL.out\" u 1:3 w linespoints lc rgb \"blue\" pt 4 t 'GMRES Full',\
\"./Times_Restarted.out\" u 1:3 w linespoints lt 1 lc rgb \"red\" pt 4 t 'GMRES Restarted'">gnuplot.gp
gnuplot gnuplot.gp
ps2pdf -dEPSCrop -dNOSAFER ./Output/Timings_GMRES.ps ./Output/Timings_GMRES.pdf
