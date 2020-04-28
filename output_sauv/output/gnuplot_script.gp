set term pngcairo font "Helvetica,18"
set output "x.png"
set xlabel 't (seconds)'
set ylabel 'x(t)'
set style fill noborder
plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title "maximal outer flowpipe", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title "maximal inner flowpipe", 'x1outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title "minimal outer flowpipe", 'x1outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,'x1inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title "minimal inner flowpipe" 
unset output
set term pngcairo font "Helvetica,18"
set output "xi.png"
set xlabel 't (seconds)'
set style fill noborder
plot 'x1outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  title "maximal outer flowpipe ", 'x1outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x1inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  title "maximal inner flowpipe ", 'x1outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 title "minimal outer flowpipe ", 'x1outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,'x1inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' title "minimal inner flowpipe ",  'x2outer.out' using 1:2 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2outer.out' using 1:3 w l lt 3 lw 2 lc rgb '#4dbeee'  notitle, 'x2inner.out' using 1:2:3 w filledcu lc rgb '#4dbeee'  notitle, 'x2outer_minimal.out' using 1:2 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle, 'x2outer_minimal.out' using 1:3 w l lt 7 lc rgb '#7e2f8e' lw 2 notitle,'x2inner_minimal.out' using 1:2:3 w filledcu lc rgb '#7e2f8e' notitle 
unset output
set term pngcairo
set output "width_ratio.png"
set xlabel 't (seconds)'
set ylabel 'min over x_i of width ratios'
plot 'width_ratio.out' w l title "width(inner-approx) / width (outer-approx)
unset output
set term pngcairo
set output "meanerror.png"
set xlabel 't (seconds)'
set ylabel 'mean over x_i of max error'
plot 'meanerror_outer.out' w l lt 1 title "outer-approx error", 'meanerror_inner.out' w l lt 1 dashtype 2 title "inner-approx error", 'meanerror_diff.out' w l lt 2 title "max distance between inner and outer-approx"
unset output
set term pngcairo
set output "relmeanerror.png"
set xlabel 't (seconds)'
set ylabel 'mean over x_i of max relative error'
plot 'relmeanerror_outer.out' w l lt 1 title "outer-approx error", 'relmeanerror_inner.out' w l lt 1 dashtype 2 title "inner-approx error", 'relmeanerror_diff.out' w l lt 2 title "max distance between inner and outer-approx"
unset output
set terminal aqua
