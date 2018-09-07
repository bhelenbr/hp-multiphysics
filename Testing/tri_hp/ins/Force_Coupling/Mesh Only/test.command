cd "$(dirname "$0")"

# This tests the motion components of the force coupling algorithm
# This test has not passed in a long time (2010)

if [ -e Results ]; then
	cd Results
else
	mkdir Results
	cd Results
fi
rm *

cp ../Inputs/* .

tri_mesh generate
mv data2_b0.grd naca.grd

rm data* rstrt*

tri_hp run

cd ..

opendiff Results Baseline