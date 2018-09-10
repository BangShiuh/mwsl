from bolos import parser, grid, solver
import numpy as np

def calculateEEDF(gas_list, gas_temperature, reduced_field, reduced_freq):
    np.seterr(divide='ignore', invalid='ignore')
    gr = grid.LinearGrid(0, 10., 500)
    boltzmann = solver.BoltzmannSolver(gr)

    with open('bolsigdb.dat') as fp:
        processes = parser.parse(fp)
    boltzmann.load_collisions(processes)

    for key in gas_list:
        boltzmann.target[key].density = gas_list[key]

    boltzmann.kT = gas_temperature * solver.KB / solver.ELECTRONVOLT
    boltzmann.EN = reduced_field
    boltzmann.FN = reduced_freq
    boltzmann.init()

    fMaxwell = boltzmann.maxwell(1.0)
    f = boltzmann.converge(fMaxwell, maxn=100, rtol=1e-5)
    # Calculate the mean energy according to the first EEDF
    mean_energy = boltzmann.mean_energy(f)
    
    # Set a new grid extending up to 15 times the mean energy.
    # Now we use a quadritic grid instead of a linear one.
    newgrid = grid.QuadraticGrid(0, 15 * mean_energy, 500)

    # Set the new grid and update the internal
    boltzmann.grid = newgrid
    boltzmann.init()

    # Calculate an EEDF in the new grid by interpolating the old one
    finterp = boltzmann.grid.interpolate(f, gr)

    # Iterate until we have a new solution
    f1 = boltzmann.converge(finterp, maxn=200, rtol=1e-5)
    Te = boltzmann.mean_energy(f1) * solver.ELECTRONVOLT / solver.KB * 2. / 3.
    if gas_temperature > Te:
        f1 = boltzmann.maxwell(gas_temperature * solver.KB / solver.ELECTRONVOLT)
        Te = gas_temperature
    mun = boltzmann.mobility(f1)
    diffn = boltzmann.diffusion(f1)
    k = boltzmann.rate(f1, "H3O+ -> H3O")

    return diffn, mun, Te, k
