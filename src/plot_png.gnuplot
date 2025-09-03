thetaxiOfile(r,norb)="./gene_data/theta_xi_O_n_orb_".norb."_J_U_".r.".dat"

rvalues="0.0 0.05 0.25"
do for [norb=2:5] {
set terminal png size 900,600 linewidth 4 font ",24"
set output "./gene_plots/plot_xi_O_norb".norb."_JU_".rvalues.".png"
# set terminal qt  linewidth 4 font ",24"
set lmargin at screen 0.16
set rmargin at screen 0.96
set tmargin at screen 0.98
set bmargin at screen 0.13
set xrange [0:0.5]
set autoscale  y
set key at 0.4,-0.05
set xlabel "{/Symbol x}" offset 0,0.5
set ylabel "O{/Sans=25({/Symbol x})}" font "Zapf Chancery, 35" offset 1,0
plot \
     for [r in rvalues] thetaxiOfile(r,norb) u 2:3 w l   t "J/U=".r
set terminal qt  linewidth 4 font ",24"
}


tildeOfile(r,norb)="./gene_data/Delta_tildeO_xi_fancyFtilde2_fancyO_n_orb_".norb."_J_U_".r.".dat"

set terminal qt  linewidth 4 font ",24"
do for [norb=2:5] {
set terminal png size 900,600 linewidth 4 font ",24"
set output "./gene_plots/plot_tildeO_Delta_norb".norb."_JU_".rvalues.".png"
set xrange [0.02:0.25]
set autoscale y
# set yrange [0:10] 
set xlabel "{/Symbol D}" offset 0,0.5
set ylabel 'ln|{~O{1\~}} |'  offset 2.0,0 font ",S"
set key at screen 0.8,0.8
plot \
     for [r in rvalues] tildeOfile(r,norb) u 1:(log(abs($2)))  w  l t "J/U=".r
set terminal qt  linewidth 4 font ",24"
}


tildeOwithDerivativesfile(r,norb)="./gene_data/Delta_tildeO_lnAbstildeO_xi_fancyFtilde2_fancyO_dtildeOdDelta_dtildeOdDeltaspline_dlnAbstildeOdDelta_dlnAbstildeOdDeltaspline_n_orb_".norb."_J_U_".r.".dat"
# we can also get Deltac from the corresponding file
get_Deltac(r,norb)=real(system("sed -n '1p' ./gene_data/Deltac_n_orb_".norb."_J_U_".r.".dat"))

# print(get_Deltac("0.0",2))


# r="0.25"
# r="0.05"
# norb="2"
do for [norb=2:5] {
do for [r in "0.0 0.05 0.25"]{
set terminal png size 900,600 linewidth 4 font ",24"
set output "./gene_plots/plot_tildeO_and_derivative_in_log_norb".norb."_JU_".r.".png"
# set border
set xrange [0.00:0.25]
set yrange [-10:8]
set xlabel "{/Symbol D}" offset 0,0.5
set ylabel 'ln|{~O{1\~}} | or dln|{~O{1\~}} |d{/Symbol D} '  offset 2.0,0 font ",S"
deltac=get_Deltac(r,norb)
plot \
     tildeOwithDerivativesfile(r,norb) u 1:3 w l t "",\
     tildeOwithDerivativesfile(r,norb) u 1:9 w l t "",\
     tildeOwithDerivativesfile(r,norb) u 1:10 w p t "",\
     sprintf("< printf '%f -10\n %f 8'",deltac,deltac) w l lt "dashed" t ""
set terminal qt  linewidth 4 font ",24"
}}




deltac(r,norb)=real(system("head -n 1 ./gene_data/Deltac_n_orb_".norb."_J_U_".r.".dat"))
deltacprecise(r,norb)=real(system("head -n 1 ./gene_data/Deltac_precise_n_orb_".norb."_J_U_".r.".dat"))
atomicO(norb,r)=0.25*norb*(r*(-norb)+r-1)
fancyF2ratio(deltapara,xi)=Power(1+Sqrt(1-(4*Power(xi,2))/Power(1-2*deltapara,2)),2)/(4.*Power(1-4*Power(xi,2),2))
# 1 theta, 2 xi, 3 O (fancyO)
thetaxiOfile(r,norb)="./gene_data/theta_xi_O_n_orb_".norb."_J_U_".r.".dat"
Power(x,y)=x**y
Sqrt(x)=sqrt(x)
maxxi(delta)=0.5-delta
gutzxi(delta)= Sqrt(1 - 4*delta)/2.0


deltacinc=0.005
do for [norb=2:5] {
do for [r in "0.0 0.05 0.25"]{
deltacval=deltacprecise(r,norb)
print(deltacval)
set terminal png size 900,600 linewidth 4 font ",24"
set output "./gene_plots/plot_fancyL_norb".norb."_JU_".r.".png"
# set terminal qt  linewidth 4 font ",24"
set lmargin at screen 0.15
set rmargin at screen 0.96
set tmargin at screen 0.98
set bmargin at screen 0.11
set xrange [0:0.25]
set autoscale y
set ytics autofreq
set ylabel "L{/Sans=25({/Symbol D},{/Symbol x})}" font "Zapf Chancery, 35" offset 1,0
set xlabel "{/Symbol x}"
deltacval1=deltacval-deltacinc
deltacval2=deltacval+deltacinc
array deltas[3]=[deltacval1,deltacval,deltacval2]
plot \
     for [i=1:3] thetaxiOfile(r,norb) u ($2):(($2>maxxi(deltas[i]))? 1/0: (-$3/atomicO(norb,r)*fancyF2ratio(deltas[i],$2)/fancyF2ratio(deltas[i],0))) w l  lc "red"  t "",\
     -1.0+0.0*x lc "black" t ""
set terminal qt  linewidth 4 font ",24"
}}



filenameresults(r,norb,tag)="./gene_data/Delta_a_b_A_U_O_E_K_Eloc_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"

set terminal png size 900,600 linewidth 4 font ",24"
set output "./gene_plots/plot_U_Delta_various_Norb_JU.png"
set xrange [0.0:0.25]
set  yrange [0.0:10.0]
set xlabel "{/Symbol D}" offset 0,1.0
set ylabel 'U/t '  offset 2.0,0 font ",S"
array norbs[3]=["2","2","5"]
array rs[3]=["0.0","0.25","0.25"]
array colors[3]=["red","black","blue"]
tag="bethe"
set key at 0.13,9
plot \
     for [i=1:3] filenameresults(rs[i],norbs[i],tag) u 1:5 w l  t "N_{orb}=".norbs[i]." J/U=".rs[i] lc rgb colors[i]
set terminal qt  linewidth 4 font ",24"
