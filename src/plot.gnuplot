# plot the data generated from example.jl


# fancyO

# 1 theta, 2 xi, 3 O (fancyO)
thetaxiOfile(r,norb)="./gene_data/theta_xi_O_n_orb_".norb."_J_U_".r.".dat"



rvalues="0.0 0.05 0.25"
do for [norb=2:5] {
set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_xi_O_norb".norb."_JU_".rvalues.".pdf"
# set terminal qt  linewidth 4 font ",24"
set lmargin at screen 0.12
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


# tildeO, i.e A^4*tildeO gives the physical O
# 1, Delta, 2 tildeO,3 xi, 4, fancyFtilde^2, 5 fancyO
tildeOfile(r,norb)="./gene_data/Delta_tildeO_xi_fancyFtilde2_fancyO_n_orb_".norb."_J_U_".r.".dat"

set terminal qt  linewidth 4 font ",24"
do for [norb=2:5] {
set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_tildeO_Delta_norb".norb."_JU_".rvalues.".pdf"
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

# we also plot the derivatives

# 1 Delta 2 tildeO 3 lnAbstildeO
# 4 xi 5 fancyFtilde2 6 fancyO
# 7  dtildeOdDelta 8 dtildeOdDeltaspline 9  dlnAbstildeOdDelta 10 dlnAbstildeOdDeltaspline
tildeOwithDerivativesfile(r,norb)="./gene_data/Delta_tildeO_lnAbstildeO_xi_fancyFtilde2_fancyO_dtildeOdDelta_dtildeOdDeltaspline_dlnAbstildeOdDelta_dlnAbstildeOdDeltaspline_n_orb_".norb."_J_U_".r.".dat"

# we can also get Deltac from the corresponding file

get_Deltac(r,norb)=real(system("sed -n '1p' ./gene_data/Deltac_n_orb_".norb."_J_U_".r.".dat"))

# print(get_Deltac("0.0",2))


# r="0.25"
# r="0.05"
# norb="2"
do for [norb=2:5] {
do for [r in "0.0 0.05 0.25"]{
set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_tildeO_and_derivative_in_log_norb".norb."_JU_".r.".pdf"
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

# landau plots

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


# r="0.25"
# deltacinc=0.002
# r="0.0"
# deltacinc=0.01
# r="0.05"
# deltacinc=0.002
# norb="2"
deltacinc=0.005
do for [norb=2:5] {
do for [r in "0.0 0.05 0.25"]{
deltacval=deltacprecise(r,norb)
print(deltacval)
set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_fancyL_norb".norb."_JU_".r.".pdf"
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


# we now explore the data for the two peak density of states

# 1 Delta 2 tildeO 3 lnAbstildeO
# 4 xi 5 fancyFtilde2 6 fancyO
# 7  dtildeOdDelta 8 dtildeOdDeltaspline
# 9  dlnAbstildeOdDelta 10 dlnAbstildeOdDeltaspline
# 11 Otp, 12 dOtpdDelta, 13 dOtpdDeltaSplin

tildeOwithDerivativesOtpfile(r,norb)="./gene_data/Delta_tildeO_lnAbstildeO_xi_fancyFtilde2_fancyO_dtildeOdDelta_dtildeOdDeltaspline_dlnAbstildeOdDelta_dlnAbstildeOdDeltaspline_Otp_dOtpdDelta_dOtpdDeltaspline_n_orb_".norb."_J_U_".r.".dat"

# r="0.0"
# r="0.05"
# r="0.25"
# norb="2"
do for [norb=2:5] {
do for [r in "0.0 0.05 0.25"]{
set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_Otp_and_dOtpdDelta_norb".norb."_JU_".r.".pdf"
set lmargin at screen 0.12
set rmargin at screen 0.96
set tmargin at screen 0.98
set bmargin at screen 0.11
set xrange [0.00:0.25]
set yrange [-5:0]
set xlabel "{/Symbol D}" offset 0,1.0
set ylabel 'O_{tp} or dO_{tp}/d{/Symbol D} '  offset 0.0,0 font ",S"
plot \
     tildeOwithDerivativesOtpfile(r,norb) u 1:11 w l t "",\
     tildeOwithDerivativesOtpfile(r,norb) u 1:12 w l t "",\
     atomicO(norb,r)+0.0*x w l t "" 
set terminal qt  linewidth 4 font ",24"
}}


# we now plot the results

# 1 Delta, 2 a, 3 b, 4 A
# 5 U, 6 O, 7 E, 8 K , 9 Eloc
filenameresults(r,norb,tag)="./gene_data/Delta_a_b_A_U_O_E_K_Eloc_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"

set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_U_Delta_various_Norb_JU.pdf"
set xrange [0.0:0.25]
set  yrange [0.0:10.0]
set xlabel "{/Symbol D}" offset 0,0.5
set ylabel 'U/t '  offset 2.0,0 font ",S"
array norbs[3]=["2","2","5"]
array rs[3]=["0.0","0.25","0.25"]
array colors[3]=["red","black","blue"]
tag="bethe"
set key at 0.13,9
plot \
     for [i=1:3] filenameresults(rs[i],norbs[i],tag) u 1:5 w l  t "N_{orb}=".norbs[i]." J/U=".rs[i] lc rgb colors[i]
set terminal qt  linewidth 4 font ",24"





get_Uc(r,norb,tag)=real(system("sed -n '1p' ./gene_data/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"))
get_Deltac(r,norb,tag)=real(system("sed -n '2p' ./gene_data/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"))
get_DeltacM(r,norb,tag)=real(system("sed -n '3p' ./gene_data/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"))
get_DeltacM1(r,norb,tag)=real(system("sed -n '4p' ./gene_data/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"))
get_DeltacI1(r,norb,tag)=real(system("sed -n '5p' ./gene_data/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"))
get_xic(r,norb,tag)=real(system("sed -n '6p' ./gene_data/Uc_Deltac_DeltacM_DeltacM1_DeltacI1_xic_lattice_".tag."_n_orb_".norb."_J_U_".r.".dat"))


# print(get_Uc("0.25","2","bethe"))
# print(get_xic("0.0","2","bethe"))


set terminal qt  linewidth 4 font ",24"
# we now add the critical points
set terminal pdfcairo size 9,6 linewidth 4 font ",24"
set output "./gene_plots/plot_U_Delta_various_Norb_JU_with_transition_line.pdf"
set xrange [0.0:0.25]
set  yrange [0.0:10.0]
set xlabel "{/Symbol D}" offset 0,0.5
set ylabel 'U/t '  offset 2.0,0 font ",S"
array norbs[3]=["2","2","5"]
array rs[3]=["0.0","0.25","0.25"]
array colors[3]=["red","black","blue"]
tag="bethe"
set key at 0.13,9
array deltacM1s[3]
array deltacI1s[3]
array deltacMs[3]
array deltacs[3]
do for [i=1:3]{deltacI1s[i]=(get_DeltacI1(rs[i],norbs[i],tag))}
do for [i=1:3]{deltacM1s[i]=(get_DeltacM1(rs[i],norbs[i],tag))}
do for [i=1:3]{deltacMs[i]=(get_DeltacM(rs[i],norbs[i],tag))}
do for [i=1:3]{deltacs[i]=(get_Deltac(rs[i],norbs[i],tag))}
plot \
     for [i=1:3] filenameresults(rs[i],norbs[i],tag) u 1:(($1>deltacs[i]||$1<deltacMs[i])?$5:1/0) w l  t "N_{orb}=".norbs[i]." J/U=".rs[i] lc rgb colors[i],\
     for [i=1:3] filenameresults(rs[i],norbs[i],tag) u 1:(($1<deltacs[i]&&$1>deltacMs[i])?$5:1/0) w l lc rgb colors[i] dashtype 2 t "",\
     for [i=1:3] sprintf("< printf '%f %f \n %f %f'", get_DeltacM1(rs[i],norbs[i],tag), get_Uc(rs[i],norbs[i],tag),get_DeltacI1(rs[i],norbs[i],tag), get_Uc(rs[i],norbs[i],tag) )  w l t ""  lc rgb colors[i] dashtype 2
set terminal qt  linewidth 4 font ",24"

# check with numeric version

set terminal qt  linewidth 4 font ",24"
# set terminal pdfcairo size 9,6 linewidth 4 font ",24"
# set output "./gene_plots/plot_U_Delta_various_Norb_JU.pdf"
set xrange [0.0:0.25]
set  yrange [0.0:10.0]
set xlabel "{/Symbol D}" offset 0,0.5
set ylabel 'U/t '  offset 2.0,0 font ",S"
array norbs[3]=["2","2","5"]
array rs[3]=["0.0","0.25","0.25"]
array colors[3]=["red","black","blue"]
tag1="bethe"
tag2="bethe_numeric"
set key at 0.13,9
plot \
     for [i=1:3] filenameresults(rs[i],norbs[i],tag1) u 1:5 w l  t "N_{orb}=".norbs[i]." J/U=".rs[i] lc rgb colors[i],\
     for [i=1:3] filenameresults(rs[i],norbs[i],tag2) u 1:5 w p lc rgb colors[i] t ""
set terminal qt  linewidth 4 font ",24"


# now, check with two peak density of states

# check two  peak
# 1 Delta 2 tildeO 3 lnAbstildeO
# 4 xi 5 fancyFtilde2 6 fancyO
# 7  dtildeOdDelta 8 dtildeOdDeltaspline
# 9  dlnAbstildeOdDelta 10 dlnAbstildeOdDeltaspline
# 11 Otp, 12 dOtpdDelta, 13 dOtpdDeltaSplin

tildeOwithDerivativesOtpfile(r,norb)="./gene_data/Delta_tildeO_lnAbstildeO_xi_fancyFtilde2_fancyO_dtildeOdDelta_dtildeOdDeltaspline_dlnAbstildeOdDelta_dlnAbstildeOdDeltaspline_Otp_dOtpdDelta_dOtpdDeltaspline_n_orb_".norb."_J_U_".r.".dat"

set terminal qt  linewidth 4 font ",24"
# set terminal pdfcairo size 9,6 linewidth 4 font ",24"
# set output "./gene_plots/plot_U_Delta_various_Norb_JU.pdf"
set xrange [0.0:0.25]
set  yrange [0.0:10.0]
set xlabel "{/Symbol D}" offset 0,0.5
set ylabel 'U/t '  offset 2.0,0 font ",S"
array norbs[3]=["2","2","5"]
array rs[3]=["0.0","0.25","0.25"]
array colors[3]=["red","black","blue"]
tag1="two_peak"
set key at 0.13,9
plot \
     for [i=1:3] filenameresults(rs[i],norbs[i],tag1) u 1:5 w l  t "N_{orb}=".norbs[i]." J/U=".rs[i] lc rgb colors[i],\
     for [i=1:3] tildeOwithDerivativesOtpfile(rs[i],norbs[i]) u 1:(-4/$12*norbs[i]) lc rgb colors[i ]    w p t ""


