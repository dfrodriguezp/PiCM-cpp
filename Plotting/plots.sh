# python3 energy_plotter.py ../../random1.json
# python3 plot_grid.py ../../random1.json --n 50
# python3 snapshots.py ../../random1.json --name phaseSpace --n 50
# python3 energy_plotter.py ../../random2.json
# python3 plot_grid.py ../../random2.json --n 50
# python3 snapshots.py ../../random2.json --name phaseSpace --n 50

# python3 snapshots.py ../../random1.json --name phaseSpace --n 50 --vlimx 3 --vlimy 3
# python3 snapshots.py ../../random2.json --name phaseSpace --n 50 --vlimx 3 --vlimy 3
# python3 snapshots.py ../../random1.json --name space --n 50
# python3 snapshots.py ../../random2.json --name space --n 50
# python3 snapshots.py ../../random1.json --name velocities --n 50 --vlimx 3 --vlimy 3 
# python3 snapshots.py ../../random2.json --name velocities --n 50 --vlimx 3 --vlimy 3 
python3 snapshots.py ../../random1.json --name directions --n 50
python3 snapshots.py ../../random2.json --name directions --n 50
