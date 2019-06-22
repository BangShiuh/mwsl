""""
A freely-propagating, premixed methane-air flat flame with applied microwave electric field.
"""

import cantera as ct
import numpy as np
import bsolver
import csv
import argparse


def parse_args():
    argparser = argparse.ArgumentParser()

    argparser.add_argument("--e",
                           help="Electric field (in V/m)",
                           type=float, default=1.0e5)
    argparser.add_argument("--h",
                           help="Enable Ohmic heating",
                           type=int, default=1)
    argparser.add_argument("--r",
                           help="update recombination",
                           type=int, default=1)
    argparser.add_argument("--t",
                           help="update elecctron Transport",
                           type=int, default=1)
    args = argparser.parse_args()
    return args


def updateElectronProperties():
    print("update electron properties")
    # Global variables
    global k
    global Te
    global mu
    global diff
    k = []
    Te = []
    mu = []
    diff = []
    # local variable
    k1 = []
    gas_list = {'N2':0,'O2':0,'H2':0,'C3H8':0,'CO2':0,'CO':0,'H2O':0}
    for n in range(f.flame.n_points):
        f.set_gas_state(n)
        gas_density = p / f.T[n] / kb
        reduced_field = np.abs(Efield/gas_density)
        reduced_freq = frequency*2*np.pi/gas_density
        for key in gas_list:
            gas_list[key] = gas.X[gas.species_index(key)]

        diffn, mun, Temp_e, rate1= bsolver.calculateEEDF(gas_list, f.T[n], reduced_field, reduced_freq)
        mu.append(mun/gas_density)#
        diff.append(diffn/gas_density)
        k1.append(rate1)
        Te.append(Temp_e)

    k = k1 + k1
    f.flame.set_electronTemperature(f.T,Te)
    if args.r:
        f.flame.set_plasmaRateCoeff(f.T, k)
    if args.t:
        f.flame.set_electronTransport(f.T, diff, mu)


# Simulation parameters
args = parse_args()
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'C3H8:0.6, O2:5, N2:18.8'  # premixed gas composition
width = 0.1  # m
loglevel = 1  # amount of diagnostic output (0 to 8)
kb = 1.38065e-23
ElectronCharge = 1.602176565e-19
Efield = args.e # [V/m]
frequency = 2.45e9 # [Hz]
multiplier = 0

# IdealGasMix object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('/home/bang/Dropbox/cantera/data/propane_ion.cti')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.IonFlame(gas, width=width)
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1, prune=0.01)
f.set_grid_min(1e-15)  #1e-10 default
f.set_min_time_step(1e-21)  #1.0e-16 default
f.max_time_step_count = 1000
f.flame.set_electricBoundaryCondition("field", 0.0)
f.radiation_enabled = False

###########################################
#                 step one                #
###########################################
f.solve(loglevel=loglevel, auto=True)

###########################################
#                 step two                #
###########################################
# update electron properties
updateElectronProperties()

# increase resolution
f.set_refine_criteria(ratio=2.5, slope=0.042, curve=0.1, prune=0.01)
# Dissociation recombination rate
if args.r:
    gas.set_multiplier(0.0, 472)
    gas.set_multiplier(0.0, 476)
    gas.set_multiplier(0.0, 478)
    #gas.set_multiplier(0.0, 482)

try:
    for i in np.linspace(0,1.0,num = 6):
        print("set plasma multiplier to  %.3e" % i)
        f.flame.set_plasmaMultiplier(i)
        f.solve(loglevel=loglevel, stage=3, enable_energy=True)
except:
    print("except")
    f.restore_steady_solution()
    f.write_csv_ND('number_density.csv', quiet=False)

print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))
SL0 = f.u[0]

###########################################
#                 step three              #
###########################################
if args.h:
    print("Begin ohmic heat")
    f.flame.set_plasmaMultiplier(1.0)
    try:
        for i in np.linspace(0,1.0,num = 15):
            print("set multiplier to  %.3e" % i)
            updateElectronProperties()
            f.flame.set_ohmicHeatingElectricField(Efield*i)
            f.solve(loglevel=loglevel, stage=3, enable_energy=True)
    except:
        print("except")
        f.restore_steady_solution()
        f.write_csv_ND('number_density.csv', quiet=False)

###########################################
#                 step four               #
###########################################
if args.h:
    for i in range(10):
        SL1 = f.u[0]
        updateElectronProperties()
        f.solve(loglevel=loglevel, stage=3, enable_energy=True)
        if (abs(f.u[0] - SL1) / f.u[0]) < 1e-3:
            print("converge")
            break

############################################
#              post simulation             #
############################################
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.u[0]))
SL1 = f.u[0]
ratio = (SL1 - SL0) / SL0
print(ratio)

updateElectronProperties()
if args.h:
    with open('flamespeed'+str(Efield), 'w') as output:
        output.write('%s' % ratio)

    with open('ohmic_heating.csv','w') as csvfile:
        writer = csv.writer(csvfile)
        for n in range(f.flame.n_points):
            f.set_gas_state(n)
            Q = Efield * Efield * mu[n] * gas.X[gas.species_index("E")] * p / kb / f.T[n] * ElectronCharge
            writer.writerow((f.grid[n],Q))

with open('chem_heat.csv','w') as csvfile2:
   writer = csv.writer(csvfile2)
   for n in range(f.flame.n_points):
       writer.writerow((f.grid[n],f.heat_release_rate[n]))

f.save('detail.xml', 'ion', 'solution with ionized gas transport')
# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('mole_fraction.csv', quiet=False)
f.write_csv_ND('number_density.csv', quiet=False)
