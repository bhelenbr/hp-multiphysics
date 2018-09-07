cd "$(dirname "$0")"

# This tests the motion components of the force coupling algorithm
# This test has not passed in a long time.
# It has never passed with petsc

HP="$HOME/Codes/tri_hp/build/Release/tri_hp_petsc"
#HP="$HOME/Codes/tri_hp/build/Release/tri_hp"

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

tri_mesh generate
mv rstrt1_b0.grd circle.grd

rm data* rstrt*

${HP} run

cd ..

opendiff Results Baseline