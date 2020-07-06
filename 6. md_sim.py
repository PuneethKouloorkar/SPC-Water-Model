from particles import *


class Simulation():

	def __init__(self, steps, input_file, box_length):
		self.steps = steps
		self.gen_input(input_file,box_length)
		self.particles = particles(input_file, 0.0001, 300, 10)
		self.dump = open('dump_file','w')
		self.energies_file = open('energy','w')
		self.trajectory_file = open('trajectory.xyz','w')
		self.energies_file.write('Timestep' + ' ' + 'Total_Energy' + ' ' + 'Temperature' '\n')
		self.dump_frequency = 10  # outputs every n steps 	

	def gen_input(self, filename, box_length):
		min_distance = 3
		coordinates = open(filename,'w')
		atom_l = int(box_length/min_distance)
		water = ([1.2361419,1.0137761,-0.0612424],[0.5104418,0.8944555,0.5514190],[1.9926927,1.1973129,0.4956931])
		masses = ([16,1,1]) 
		x_vec = np.transpose([1,0,0])
		y_vec = np.transpose([0,1,0])
		z_vec = np.transpose([0,0,1])
		atom_count = 0 
		molecule_count = 1

		for x in range(atom_l):
			for y in range(atom_l):
				for z in range(atom_l):

					opt_water = water 
					opt_water = opt_water + x*min_distance*x_vec + y*min_distance*y_vec + z*min_distance*z_vec

					for i in range(len(water)):
						coordinates.write('%s %s %s %s %s %s \n' % (atom_count, masses[i], molecule_count, opt_water[i][0],opt_water[i][1],opt_water[i][2]))
						atom_count += 1
						if atom_count % len(water) == 0:
							molecule_count += 1

		coordinates.close()


	def close_files(self):

		self.dump.close()
		self.energies_file.close()	
		self.trajectory_file.close()	

	def run(self):

		self.particles.centre_velocity()
		self.particles.calc_intforces()
		self.particles.calc_extforces()
		

		for t in range(self.steps):
			self.update(t) 

			#print("step " + str(t))
		self.close_files()

	def update(self, timestep):

		if timestep % self.dump_frequency == 0: 
			self.total_energy, self.temperature = self.particles.calculate_hamiltonian() 
			self.particles.write_output(timestep,self.dump,self.total_energy,self.temperature)
			

		self.particles.write_energies(timestep,self.energies_file)
		self.particles.write_trajectory(timestep,self.trajectory_file)
		
		self.particles.update_velocities() # intermediate velocity, t+1/2dt vel-verlet

		self.particles.update_positions()
		self.particles.PBC()
		
		self.particles.calc_intforces()
		self.particles.calc_extforces()
		self.particles.update_acceleration()
		self.particles.update_velocities() # computed velocity, with update acceleration at t + dt
		self.particles.calculate_hamiltonian()