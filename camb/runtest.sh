#!/bin/env bash

make clean
make


./camb params_GP.ini
./camb params_thetabin.ini
./camb params_smoothbin.ini

python coupling_plot.py 
python matpower_plot.py 
python solutions_plot.py  
python spectra_plot.py
