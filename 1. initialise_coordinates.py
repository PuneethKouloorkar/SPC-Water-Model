import numpy as np 


def initialise_positions(box_length, min_distance, pos_file):
	coordinates = open(pos_file,'w')
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
				print(water)
				opt_water = opt_water + x*min_distance*x_vec + y*min_distance*y_vec + z*min_distance*z_vec
				# print(min_distance*x_vec,opt_water)
				for i in range(len(water)):
					coordinates.write('%s %s %s %s %s %s \n' % (atom_count, masses[i], molecule_count, opt_water[i][0],opt_water[i][1],opt_water[i][2]))
					atom_count += 1
					if atom_count % len(water) == 0:
						molecule_count += 1
		
			

	coordinates.close()

initialise_positions(12,3,'positions')