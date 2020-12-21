#!/bin/sh
grep -A 0 "Currents" mos_tran0.5_19_eafe1.txt > temp_mos.txt
#!/bin/awk
#awk -F '[ ,]' '/anode/{print $10}' temp_semi.txt > current_eafe1_reverse.txt
awk -F '[ ,]' '/drain/{print $12}' temp_mos.txt > current_mos_tran0.5_19_eafe1.txt
#awk -F '[ ,]' '/emitter/{print $12}' temp_pnp.txt > current_pnp_eafe1_CB0.5.txt
#awk -F '[ ,]' '/emitter|base|collector/{print $12, $15, $18}' temp_pnp.txt > BJT1/current_pnp_eafe1_BE0.02.txt

