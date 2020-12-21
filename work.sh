make
bsub -J mos_tran0.5_19_eafe1_NOT12 -q batch -R "span[ptile=36]" -n 144 -o "mos_tran0.5_19_eafe1_NOT12.txt" "mpiexec ./semiconductor"
#bsub -J npn_eafe_CB3 -q batch -R "span[ptile=18]" -n 72 -o "npn_eafe_CB3.txt" "mpiexec ./semiconductor"
#bsub -J npn_eafe_CE2.5 -q batch -R "span[ptile=18]" -n 72 -o "npn_eafe_CE2.5.txt" "mpiexec ./semiconductor"
#bsub -J pnp_eafe_CE1 -q batch -R "span[ptile=36]" -n 144 -o "pnp_eafe_CE1_test.txt" "mpiexec ./semiconductor"
#bsub -J pnp_eafe_BE0.02 -q batch -R "span[ptile=18]" -n 72 -o "pnp_eafe_BE0.02.txt" "mpiexec ./semiconductor"
#bsub -J pnp_eafe_CB0.5 -q batch -R "span[ptile=18]" -n 72 -o "pnp_eafe_CB0.5_big.txt" "mpiexec ./semiconductor"
#bsub -J mos3.4_18_eafe1 -q batch -R "span[ptile=36]" -n 144 -o "mos3.4_18_eafe1_vtk.txt" "mpiexec ./semiconductor"
#bsub -J semi_eafe1 -q batch -R "span[ptile=36]" -n 144 -o "semi_eafe1_vtk.txt" "mpiexec ./semiconductor"
