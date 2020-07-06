import numpy as np 
import math   

DIMS = 3

KB = 0.0083  # kcal /mol K 

class particles():
	def __init__(self, input_file, step_length, set_temperature, box_length):
	
		self.step_length = step_length 
		self.set_temperature = set_temperature # for eventual thermostatting
		self.v_int = 0
		self.v_ext = 0
		self.box = box_length
		self.mass = np.genfromtxt(input_file, usecols=[1])
		self.bond_data = np.genfromtxt(input_file, usecols=[2])  # cols 0-2 contain number no, type and bonds 
		self.coordinates = np.genfromtxt(input_file, usecols=[3,4,5])
		self.n_particles = len(self.bond_data)
		self.molecules = int(self.n_particles/3) 
		self.variance = np.sqrt(KB*set_temperature/self.mass)
		self.maxdist = np.random.normal(0,1,(DIMS,self.n_particles))
		self.velocities = np.transpose(self.variance*self.maxdist) #  gaussian velocities with sigma^2 = init_vel
		self.acceleration = np.zeros((self.n_particles,DIMS))  # initialises N x 3 array of zeroes
		self.forces = np.zeros((self.n_particles,DIMS)) # initialises N x 3 array of zeroes
		self.forceOH = 1054.2     # force constant kcal mol-1 A^-2
		self.forcetheta = 75.9  # force constant kcal mol-1 radian-2
		self.bondOH = 1.0    # r0 in angstroms  
		self.thetaOH = 1.91  # theta0 in radians theta0 = 109,47 deg
		

	def centre_velocity(self):
		mass_vel = np.dot(self.mass,np.linalg.norm(self.velocities,axis=1))
		v_cm = np.divide(mass_vel,np.sum(self.mass)) # calculates velocity centre of mass 
		self.velocities = self.velocities - v_cm

	def update_positions(self):
		#print(str(self.coordinates))
		for i in range(self.n_particles):
			for j in range(DIMS):
				self.coordinates[i][j] = self.coordinates[i][j] + self.step_length*self.velocities[i][j] + ((self.step_length**2)/2)*self.acceleration[i][j]
				
	def PBC(self):
		for i in range(self.molecules):
			for j in range(DIMS):
				if self.coordinates[i*3][j] < -0.5:	
					#print(str(self.coordinates[i][j])+' first under 0')
					self.coordinates[i*3][j] += (self.box -0.5)
					#print(str(self.coordinates[i][j]))
					self.coordinates[i*3+1][j] += (self.box - 0.5)
					self.coordinates[i*3+2][j] += (self.box -0.5)

				if self.coordinates[i*3][j] > self.box + 0.5:
					#print(str(self.coordinates[i][j])+' first over box')
					self.coordinates[i*3][j] -= (self.box + 0.5)
					#print(str(self.coordinates[i][j]))
					self.coordinates[i*3+1][j] -= (self.box +0.5)
					self.coordinates[i*3+2][j] -= (self.box +0.5)

	def update_velocities(self):
		for i in range(self.n_particles):
			for j in range(DIMS):
				self.velocities[i][j] = 0.5*self.step_length*self.acceleration[i][j] + self.velocities[i][j]

	
	def velocity_squared(self,velocities):
		velocity_sq3D = np.multiply(velocities,velocities)
		velocity_sq  = np.sum(velocity_sq3D,axis=1)

		return velocity_sq

	def update_acceleration(self):
		for i in range(self.n_particles):
			for j in range(DIMS):
				self.acceleration[i][j] = self.forces[i][j] / self.mass[i]  # force currently 1, gets mass from column 1

	def calculate_ekin(self):
		velocity_sq = self.velocity_squared(self.velocities)/3   # need average over three axes

		energy_kinetic = 0.5*np.dot(self.mass,velocity_sq)   # 1/2 sum mivi^2 over x, y ,z
		
		return energy_kinetic 

	def calculate_temperature(self):
		energy_kinetic = self.calculate_ekin()
		temperature = (2*energy_kinetic)/(3*self.n_particles*KB) #quipartition theorem, 3N (3N-3 minus vibrations), minus velocities -->3/2 

		return temperature 

	def difference_vector(self,p1,p2):
		difference = p2-p1

		return difference

	def calc_magnitude(self,vec):
		if np.linalg.norm(vec) < 0.5*self.box:
			norm = np.linalg.norm(vec)
			#print(str(norm) + '  1.')
		else:
			norm = np.linalg.norm(vec)-self.box
			#print(str(norm) + '  2.')
		
		return np.abs(norm)
		
	def calculate_hamiltonian(self):
		energy_kinetic, temperature = self.calculate_ekin(), self.calculate_temperature()
		hamiltonian = energy_kinetic + self.v_int + self.v_ext 

		return hamiltonian, temperature

	def write_output(self,step,file,hamiltonian,temperature):
		file_header = 'Timestep' + ' ' + 'x' + ' ' + 'y' + ' ' + 'z' + ' ' + 'Total_Energy' + ' ' + 'Temperature' '\n' 
		file.write(file_header)	

		for atom in range(self.n_particles):
			file.write('%s %s %s %s %s %s %s \n' % (atom, step, self.coordinates[atom][0], self.coordinates[atom][1], self.coordinates[atom][2], hamiltonian, temperature))
			# file.write(%s + '' %s + '' + %s + '' + %s + '' + %s + '' % (step, self.coordinates[i])'\n')

	def write_trajectory(self,step,file):
		atom_dict = {16: 'O', 1: 'H'}
		file_header = str(self.n_particles) + '\n' + 'Timestep: %s \n' % step  
		file.write(file_header)
		
		for atom in range(self.n_particles):
			file.write('%s %s %s %s \n' % (atom_dict.get(self.mass[atom]), self.coordinates[atom][0], self.coordinates[atom][1], self.coordinates[atom][2]))
		
	def write_energies(self,step,file):
		hamiltonian, temperature = self.calculate_hamiltonian()
		file.write('%s %s %s \n' % (step, hamiltonian,temperature))


	def calc_intforces(self):
		self.v_int = 0 
		self.forces = np.zeros((self.n_particles,DIMS))

		for i in range(self.molecules):

			OH1 = self.difference_vector(self.coordinates[(3*i)],self.coordinates[(3*i+1)]) # vector O-H1
			#print(str(OH1) + ' OH1')
			OH2 = self.difference_vector(self.coordinates[(3*i)],self.coordinates[(3*i+2)]) # vector O-H2
			#print(str(OH2) + ' OH2')
			bond1 = self.calc_magnitude(OH1)
			#bond1 = np.linalg.norm(OH1)
			#bond2 = np.linalg.norm(OH2)

			#print(str(bond1)+ ' bond1')
			bond2 = self.calc_magnitude(OH2)
			#print(str(bond2)+ ' bond 2')
		
			theta = np.dot(OH1,OH2)/(bond1*bond2)  # definition of dot product 

			harmonic_r2 = 0.5*self.forceOH*((bond1 - self.bondOH)**2)
			harmonic_r1 = 0.5*self.forceOH*((bond2 - self.bondOH)**2)
			harmonic_theta = 0.5*self.forcetheta*(theta - self.thetaOH)**2

			self.v_int +=  harmonic_theta + harmonic_r1 + harmonic_r2
			#self.v_int +=  harmonic_r2 + harmonic_r1
			#self.v_int +=  harmonic_theta

			unit_OH1 = OH1/bond1
			unit_OH2 = OH2/bond2

			# force(r) harmonic potential  dV/dr
			forcer_h1 = -self.forceOH*(bond1-self.bondOH)*unit_OH1   # force on H1
			forcer_h2 = -self.forceOH*(bond2-self.bondOH)*unit_OH2   # force on H2 
			forcer_O = -(forcer_h1 + forcer_h2) # force on O  

			# force(theta) harmonic potential dV/dtheta
			forcet_h1 = -self.forcetheta*(theta-self.thetaOH)*(((OH2)/(bond1*bond2))+((OH2-OH1)/bond2)*unit_OH1*(1/np.sqrt(bond1)))
			#print(forcet_h1)
			forcet_h2 = -self.forcetheta*(theta-self.thetaOH)*(((OH1)/(bond1*bond2))+((OH1-OH2)/bond1)*unit_OH2*(1/np.sqrt(bond2)))
			#print(forcet_h2)
			forcet_O = -(forcet_h1+forcet_h2)
			# print(forcer_h1,forcet_h2,forcet_O)

			# update force array 
			self.forces[3*i] =  forcer_O + forcet_O
			self.forces[3*i+1] = forcer_h1 + forcet_h1
			self.forces[3*i+2] = forcer_h2 + forcet_h2

		#print(str(self.v_int) + ' internal')

	def calc_extforces(self):
		self.v_ext = 0
		r_dist = np.zeros((self.molecules-2,DIMS))
		r_norm = np.zeros(len(r_dist)-1)
		E_LJ = 0.0
		E_LJ_F = np.zeros((self.molecules,DIMS))

		for i in range(self.molecules-2):
			for j in range(i+1,self.molecules-1,1):
				r_dist[i] = self.difference_vector(self.coordinates[i*3],self.coordinates[j*3]) #self.coordinates(x*3) is Oxygen

		#print(r_dist)		

		for i in range(len(r_dist)-1):
			r_norm[i] = np.linalg.norm(r_dist[i])
		#print(r_norm)
		
		epsilon = 0.1554253					#[kcal/mol]
		sigma = 3.165492					#Angstr.
		
		A = 4*epsilon*np.power(sigma,12)
		B = 4*epsilon*np.power(sigma,6)
		C = 332.0636						#kcal*Angstr/mol

		c1 = 0.41
		c2 = -0.82
		
		R_cutoff = 2.5*sigma
		phi_cutoff = 4*epsilon*(np.power(R_cutoff,-12)- np.power(R_cutoff,-6))

		#calc energy and forces external
		#forces are derived from Lennard Jones Pot.
		#forces on both atoms, that means on second one the vector is negative
		for i in range(len(r_norm)-2):
			for j in range(i+1,len(r_norm)-1,1):
				if r_norm[i] <= R_cutoff:
					E_LJ = ((A*np.power(r_norm[i],-12))-(B*np.power(r_norm[i],-6)) - phi_cutoff)
					E_LJ_F[i] = ((-12*A*np.power(r_norm[i],-13)+6*B*np.power(r_norm[i],-7)) - C*(c2**2*np.power(r_norm[i],-2)))*(r_dist[i]/(np.sqrt(r_norm[i])*r_norm[i]))
					
					self.forces[3*i] += E_LJ_F[i]
					#self.forces[3*i+1] += E_LJ_F[i]
					#self.forces[3*i+2] += E_LJ_F[i]
					self.forces[3*j] -= E_LJ_F[i]
					#self.forces[3*j+1] -= E_LJ_F[i]
					#self.forces[3*j+2] -= E_LJ_F[i]
				
					self.v_ext += E_LJ + C*(c2**2/r_norm[i])
				else:
					self.forces[3*i] += -C*(c2**2*np.power(r_norm[i],-2))*(r_dist[i]/(np.sqrt(r_norm[i])*r_norm[i]))
					self.forces[3*j] -= -C*(c2**2*np.power(r_norm[i],-2))*(r_dist[i]/(np.sqrt(r_norm[i])*r_norm[i]))
					self.v_ext += C*(c2**2/np.abs(r_norm[i]))

		#print(str(self.v_ext) + ' external')
		return self.v_ext