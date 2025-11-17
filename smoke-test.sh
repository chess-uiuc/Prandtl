#!/bin/bash
# This demonstrates a smoke test for Prandtl alone.
# Ideally, if the problem is coupled say HEGEL/CHyPS, or HEGEL/Flux
# Then please describe the steps needed to prepare the data to run
# It is fine to pre-prepare the data, and add it to your PR.
# The goal is that one can read the script and piece together how
# to run a case.
set -e # says exit with nonzero exit code on error
mkdir -p smoke_test
cd smoke_test
ln -s ../Prandtl .
cp -r ../TestCases/Euler/1D/AcousticWave .
cd AcousticWave
# Now EDIT config.json (edit the file paths):
#         "output_file_path" : "../../TestCases/Euler/1D/AcousticWave",
#      -->"output_file_path" : "./"
#         "mesh_file": "../../TestCases/Euler/1D/AcousticWave/AcousticWave.mesh",
#      -->"mesh_file": "./AcousticWave.mesh"
# Save as config_smoke.json
mpiexec -n 2 ../Prandtl -c ./config_smoke.json
# The run should produce output files:
# Paraview/Paraview.pvd
# Paraview/Cycle000<NNN>
# This test simply checks for the existence of the output files
if [[ -f Paraview/Paraview.pvd && -d Paraview/Cycle000000 ]]; then
    echo "Smoke test passed.  Output files exist."
    exit 0
else
    echo  "Smoke test failed: Did not produce output files."
    exit 1
fi
