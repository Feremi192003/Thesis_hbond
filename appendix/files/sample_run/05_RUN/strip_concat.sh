if [ $1 == "strip" ]; then
    cat << EOF > strip.in
parm $2
trajin $3

strip :WAT
strip :Na+

trajout $4
EOF
    cpptraj -i strip.in > log/cpptraj_strip_${5}.log
elif [ $1 == "concat" ]; then
    cat << EOF > concat.in
parm $2
trajin $3
trajin $4

center
autoimage

trajout $5
EOF
    cpptraj -i concat.in
else
    echo "Error"
fi
