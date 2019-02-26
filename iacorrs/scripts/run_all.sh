#!/bin/bash

echo "Generating blind parameter values"
python3.6 genblind.py

python3.6 -m pdb get_wgg.py --config settings.yaml --redshift 0.062 --process wgg wgp wpp --ind 0
python3.6 get_wgg.py --config settings.yaml --redshift 0.300 --process wgg wgp wpp --ind 2
python3.6 get_wgg.py --config settings.yaml --redshift 0.625 --process wgg wgp wpp --ind 4
python3.6 get_wgg.py --config settings.yaml --redshift 1.000 --process wgg wgp wpp --ind 6

python3.6 bind.py