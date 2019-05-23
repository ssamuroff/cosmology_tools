#!/bin/bash
#SBATCH --job-name=gen_blinded_data
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10000
#SBATCH -t 8:00:00
#SBATCH -o /home/ssamurof/logs/blinddata.log


echo "Generating blind parameter values"
python genblind.py

python get_wgg.py --config settings.yaml --redshift 0.062 --process wgg wgp wpp --ind 0
python get_wgg.py --config settings.yaml --redshift 0.300 --process wgg wgp wpp --ind 2
python get_wgg.py --config settings.yaml --redshift 0.625 --process wgg wgp wpp --ind 4
python get_wgg.py --config settings.yaml --redshift 1.000 --process wgg wgp wpp --ind 6

python bind.py