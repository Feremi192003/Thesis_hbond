restraint=(100 50 25 10 5 1)

for i in ${!restraint[*]}; do
	cat << EOF > equil${i}.in
2ns equilibration with restraint ${restraint[$i]}
&cntrl
  imin=0,
  ntx=5,
  irest=1,
  ntpr=2000,
  ntwx=2000,
  ntf=2,
  ntc=2,
  ntb=2,
  nstlim=2000000,
  temp0=300,
  ntt=3,
  gamma_ln=5.0,
  ntp=1,
  taup=2.0,
  dt=0.001,
  restraint_wt=${restraint[$i]},
  restraintmask=':1-277'
&end
&wt
  type='END',
&end
EOF
done

