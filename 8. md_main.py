import numpy as np 

from md_sim import *

import matplotlib.pyplot as plt

# notes on input file format : col 0 = molecule number, col 1 = mass, col 2= molecule number, col 3-5 = positions 

def main():

<<<<<<< HEAD
			if timestep % dump_frequency == 0: 
				total_energy, temperature = simulation.calculate_hamiltonian() 
				simulation.write_output(timestep,dump,total_energy, temperature)
				
			
			simulation.write_energies(timestep,energies_file)
			simulation.write_trajectory(timestep,trajectory_file)
			
			simulation.update_positions()
			simulation.update_velocities() # intermediate velocity, t+1/2dt vel-verlet
			simulation.calc_intforces()
			simulation.update_acceleration()
			simulation.update_velocities() # computed velocity, with update acceleration at t + dt
			#simulation.leap_frog()
			simulation.calculate_hamiltonian()
			simulation.calc_extforces()
 
		
	dump.close()
	energies_file.close()	
	trajectory_file.close()	
=======
	sim = Simulation(40000, 'positions', 10)	
	sim.run()
	
>>>>>>> int_forces

if __name__ == '__main__':	
	main()