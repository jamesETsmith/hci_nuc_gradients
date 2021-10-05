#!/bin/bash

JAVA_PATH=/usr/bin/java
JMOL_PATH=/mnt/home/jsmith/apps/jmol-14.31.35/Jmol.jar

mkdir -p 1_A/40e_40o 1_B/40e_40o 1_C/40e_40o 3_A/40e_40o 3_B/40e_40o 3_C/40e_40o
mkdir -p final/1_A/40e_40o final/1_B/40e_40o final/1_C/40e_40o final/3_A/40e_40o final/3_B/40e_40o final/3_C/40e_40o

python3 jmol_write_mos.py ../../dz/1_A/40e_40o/_molden/1_A_DZ_40e_40o_eps1=7.5e-05_initial.molden $JAVA_PATH $JMOL_PATH 1_A/40e_40o/
python3 jmol_write_mos.py ../../dz/1_B/40e_40o/_molden/1_B_DZ_40e_40o_eps1=7.5e-05_initial.molden $JAVA_PATH $JMOL_PATH 1_B/40e_40o/
python3 jmol_write_mos.py ../../dz/1_C/40e_40o/_molden/1_C_DZ_40e_40o_eps1=7.5e-05_initial.molden $JAVA_PATH $JMOL_PATH 1_C/40e_40o/
python3 jmol_write_mos.py ../../dz/3_A/40e_40o/_molden/3_A_DZ_40e_40o_eps1=7.5e-05_initial.molden $JAVA_PATH $JMOL_PATH 3_A/40e_40o/
python3 jmol_write_mos.py ../../dz/3_B/40e_40o/_molden/3_B_DZ_40e_40o_eps1=7.5e-05_initial.molden $JAVA_PATH $JMOL_PATH 3_B/40e_40o/
python3 jmol_write_mos.py ../../dz/3_C/40e_40o/_molden/3_C_DZ_40e_40o_eps1=7.5e-05_initial.molden $JAVA_PATH $JMOL_PATH 3_C/40e_40o/

# python3 jmol_write_mos.py ../../dz/1_A/40e_40o/_molden/1_A_DZ_40e_40o_eps1=7.5e-05_final.molden $JAVA_PATH $JMOL_PATH final/1_A/40e_40o
# python3 jmol_write_mos.py ../../dz/1_B/40e_40o/_molden/1_B_DZ_40e_40o_eps1=7.5e-05_final.molden $JAVA_PATH $JMOL_PATH final/1_B/40e_40o
# python3 jmol_write_mos.py ../../dz/1_C/40e_40o/_molden/1_C_DZ_40e_40o_eps1=7.5e-05_final.molden $JAVA_PATH $JMOL_PATH final/1_C/40e_40o
# python3 jmol_write_mos.py ../../dz/3_A/40e_40o/_molden/3_A_DZ_40e_40o_eps1=7.5e-05_final.molden $JAVA_PATH $JMOL_PATH final/3_A/40e_40o
# python3 jmol_write_mos.py ../../dz/3_B/40e_40o/_molden/3_B_DZ_40e_40o_eps1=7.5e-05_final.molden $JAVA_PATH $JMOL_PATH final/3_B/40e_40o
# python3 jmol_write_mos.py ../../dz/3_C/40e_40o/_molden/3_C_DZ_40e_40o_eps1=7.5e-05_final.molden $JAVA_PATH $JMOL_PATH final/3_C/40e_40o