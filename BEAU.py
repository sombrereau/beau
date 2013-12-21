class depletion_parms:

	# This object holds all the parameters for depletion analysis

	def __init__(self):

		# depletion parameters for fuel cycle analysis
		
		self.runTRANS           = 'yes'
		self.nBU                =     5
		self.nBU1               =     3
		self.dBU1               =  1000
		self.dBU2               = 10000
		self.library            = 'pwru50'
		self.mode               = 'source'
		self.method             = 'beginning'
		self.start_ocf          = 0
		self.format             = 'no'
		self.transport_module   = 'MCNP5'
		self.fission_power_only = 'no'

		# depletion parameters for fuel cycle perturbation

                self.burnup_err     = 0.005
		self.burnup_perturb = .2		

                self.power_err     = 0.005
		self.power_perturb = .2

                self.mat_err    = 0.01
                self.keff_err   = 0.001

		self.import_nuclides = {}
		# this is a dictionary of ORIGEN format lists of nuclides indexed by cell
		# self.important_nuclides = {'13100':['922350','922380']}

		self.nMat = 5
		# this is the minimum number of times the fuel cycle must be modeled before BEAU starts checking for convergence
			
		self.mat_weighting = .5
		# this is the weighting factor of the N+1 composition of the Nth composition
		# mat (N+1) = mat(N)*(self.mat_weighting) + mat_estimate(N+1)*(1 - self.mat_weighting)

class equilibrium_cycle:
	
	# This object holds all the information for the fuel cycle of pebblebed reactor for a specific model given the core's discitization

	def __init__(self):
		
		import time

		self.depl_parms = depletion_parms()
		self.target_burnups = {}	
		# this defines the target burnups for each circuit progression
		# self.target_burnups = {<circ_prog_key>:target burnup, ... }

		self.target_power = 0 	
		# this defines the target power for system
		self.target_k = 1.
		self.AUSCE = AUSCE()
		self.initial_input = open('inp1').read()
		self.initial_time  = time.ctime()

	def populate(self, MPO=None, cycle_length=None, neutron_source=None,):
		
		import mocup

		self.mpo  = MPO		
		self.mpos = self.mpo.split()
		self.circ_progs = {}
		self.cycle_length = cycle_length
		self.neutron_source = neutron_source

		for circ_prog_key in self.fuel_cycle.keys():
			PF = {}
			mpo1 = mocup.mpo()
			for circ_key in self.fuel_cycle[circ_prog_key].keys():
				for prog_key in self.fuel_cycle[circ_prog_key][circ_key].keys():
					PF[prog_key] = self.pebble_flow[prog_key]
					for cell in self.fuel_cycle[circ_prog_key][circ_key][prog_key]:
						mpo1 = self.mpos[cell] + mpo1 
						
			circ_prog = circuit_progression()
			circ_prog.makeup = self.makeup[circ_prog_key]
			circ_prog.populate(title=circ_prog_key,MPO=mpo1, cycle_length=self.cycle_length[circ_prog_key],fuel_cycle=self.fuel_cycle[circ_prog_key],circuits=self.circuits[circ_prog_key], depl_parms=self.depl_parms, neutron_source=self.neutron_source,burnup=self.target_burnups[circ_prog_key], pebble_flow = PF)
			self.circ_progs[circ_prog_key] = circ_prog

	def power(self):

		total_power = 0 
		for circ_prog_key in self.fuel_cycle.keys():
			pow = self.circ_progs[circ_prog_key].power()
			total_power += pow
		return total_power		

	def search(self):

		for circ_prog_key in self.fuel_cycle.keys():
			self.circ_progs[circ_prog_key].search_cycle_length()

	def deplete(self):
			
		for circ_prog_key in self.fuel_cycle.keys():
			self.circ_progs[circ_prog_key].deplete()

	def burnup(self):
		
		self.ave_BU = {}
		for circ_prog_key in self.fuel_cycle.keys():
			self.ave_BU[circ_prog_key] = self.circ_progs[circ_prog_key].burnup()
		return self.ave_BU
		
	def compBRO(self,iteration=None):
		
		#This method generates a new input equilibrium input deck based on the equilibrium fuel cycle just modeled
		
		import mocup
		import os
		import re
		
		# ensure depletion has been modeled
		
		for circ_prog_key in self.fuel_cycle.keys():
			for circ_key in self.fuel_cycle[circ_prog_key].keys():
				for prog_key in self.fuel_cycle[circ_prog_key][circ_key].keys():
					for cell in self.fuel_cycle[circ_prog_key][circ_key][prog_key]:
						if len(self.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].mpos[0].ORIGEN_power.keys()) == 0:
							print 'cell: %s in progression: %s in circuit: %s has not yet been depleted' % (cell, prog_key, circ_key)
		
		input = open('inp1').read()
		
		# progress through depletion cells
		
		k = input.find('\n\n')
		cell_strings = input[:k]
		
		cell_tokens = re.compile('\n\s{0,4}\d{1,5}').findall(cell_strings)
		
		self.new_mats = {}
		
		for circ_prog_key in self.fuel_cycle.keys():
			for circ_key in self.fuel_cycle[circ_prog_key].keys():
				for prog_key in self.fuel_cycle[circ_prog_key][circ_key].keys():
					for cell in self.fuel_cycle[circ_prog_key][circ_key][prog_key]:
						m = int(self.mpos[cell].mat_ID[cell])
			
						# find material card
						input_mat = mocup.material()
						if self.depl_parms.transport_module[0] in 'mM':
							input_mat.import_mcf('inp1',m)
						elif self.depl_parms.transport_module[0] in 'sS':
							input_mat.import_scf('inp1',m)
						token = input_mat.mat_card
			
						eq_mat = self.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].time_ave_mat()[cell]
						# normalize to at/bn-cm
						eq_mat = eq_mat*(.60221415/self.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].mpos[0].volume[cell])

						if iteration > 1.:
							# include previous estimate of eq mat in new estimate of eq mat
							old = mocup.material()
							eq_ocf_loc = 'moi_files/moi.%s.eq.pch'
							old.import_ocf(eq_ocf_loc)
							
							new = eq_mat*(1.-self.mat_weighting) + old*self.mat_weighting
							new.make_ocf(eq_ocf_loc,lib=self.depl_parms.library)

							for nucl in input_mat.comp.keys():
								# this loop over rides concentrations of nuclides that are depleted
								if nucl in new.comp.keys():
									input_mat.comp[nucl] = new.comp[nucl]
						else:
							for nucl in input_mat.comp.keys():
								# this loop over rides concentrations of nuclides that are depleted
								if nucl in eq_mat.comp.keys():
									input_mat.comp[nucl] = eq_mat.comp[nucl]

						if cell in self.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].mpos[0].temperature.keys():
							if self.depl_parms.transport_module[0] in 'mM':
								mat_card = input_mat.mcf(m,tmp=self.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].mpos[0].temperature[cell])
						elif new.tmp != 0:
							if self.depl_parms.transport_module[0] in 'mM':
								mat_card = input_mat.mcf(m)
						else:
							if self.depl_parms.transport_module[0] in 'mM':
								mat_card = input_mat.mcf_mcnp_nuclides(m)
						if self.depl_parms.transport_module[0] in 'sS':
							mat_card = input_mat.scf(m)

						input = input.replace(token,mat_card)
						self.new_mats[cell] = input_mat		

		if self.depl_parms.transport_module[0] in 'mM':
			for token in cell_tokens:
				cell = token.replace('\n','').replace(' ','')
				if cell in self.new_mats.keys():
					k = cell_strings.find(token) + len(token)
					l = cell_strings[k:].find(self.mpos[cell].mat_ID[cell]) + len(self.mpos[cell].mat_ID[cell])
					n = cell_strings[k+l:].find(cell_strings[k+l:].split()[0])+len(cell_strings[k+l:].split()[0])
					
					token1 = token + cell_strings[k:k+l+n]
					new_token1 = '\n%s %s %1.5E' % (cell, self.mpos[cell].mat_ID[cell], self.new_mats[cell].moles())
					input = input.replace(token1,new_token1)
		
		open('inp1','w').write(input)
				
	def search_source(self, ):

		import math
		from os import system
		
		system('\\rm scratch source.o source.csv')
		system('touch scratch/source.o')
		system('touch scratch/source.csv')
		csv  = 'neutron source search, target power: %1.5E\niterations (#),source (n/s),power (MWt)' % (self.target_power)
		open('scratch/source.csv','w').write(csv)
		
		# search for neutron source to impose a target power
		
		neutron_source = self.target_power/(200.*1.6022e-19)*2.5
		
		powers = {(1.-self.depl_parms.power_perturb)*neutron_source:0, (1.+self.depl_parms.power_perturb)*neutron_source:0}	

		l = 0 

		for NS in powers.keys():
			l += 1
			self.populate(MPO=self.mpo, cycle_length=self.cycle_length, neutron_source=NS)
			csv = '\n%d,%1.5E,' % (l,NS)
			open('scratch/source.csv','a').write(csv)
			# deplete all the circuit progessions in the equilibrium fuel cycle
			self.search()
			csv = '%1.5E' % self.power()
			print csv
			open('scratch/source.csv','a').write(csv)
			powers[NS] = self.power()
			system('cat scratch/cycle_length.o >> scratch/source.o')

		min = 0 
		max = 0 
		
		while math.fabs(self.power()/self.target_power-1.) > self.depl_parms.power_err:
			
			l += 1
			
			m,b = least_squares(powers)
			
			NS = (self.target_power - b)/m

			# ensure neutron source is positive

			if NS < .1*sorted(powers.keys())[0]:
				NS = .1*sorted(powers.keys())[0]
			
			elif NS > 10*sorted(powers.keys())[-1]:
				NS = 10*sorted(powers.keys())[-1]
			
			elif (min > 6) or (max > 6):
				NS = 0.5*(sorted(powers.keys())[0] + sorted(powers.keys())[-1])
			
			self.populate(MPO=self.mpo, cycle_length=self.cycle_length, neutron_source=NS)
			csv = '\n%d,%1.5E,' % (l,NS)
			open('scratch/source.csv','a').write(csv)
			self.search()
			csv = '%1.5E' % self.power()
			open('scratch/source.csv','a').write(csv)
			powers[NS] = self.power()
			system('cat scratch/cycle_length.o >> scratch/source.o')
			
			if math.fabs(self.power()/self.target_power-1) < self.depl_parms.power_err:
				# neutron source has converged!!
				pass
			else:
				# ensure best points for interpolation / extrapolation
				errs = {}
				for NS in powers.keys():
					errs[math.fabs(powers[NS] - self.target_power)] = NS
				del powers[errs[sorted(errs.keys())[-1]]]
				
				if self.power() > self.target_power:
					min  = 0 
					max += 1
				else:
					max  = 0 
					min += 1
		
		system('touch scratch/source.csv')
		system('cat scratch/source.csv >> source.o')

class circuit_progression:

	# This object holds all the information for the fuel cycle of pebblebed reactor for a specific model given the core's discitization for a given feed fuek

	def __init__(self):

		import mocup

		self.passes = {}
		self.makeup = mocup.material()
		self.target_burnup = 0 

	def populate(self,title=None, MPO=None, cycle_length=None, neutron_source=None, fuel_cycle=None, circuits=None, depl_parms=None, burnup=None, pebble_flow=None):
	
		import mocup

		self.mpo = MPO
		self.title =title
		self.cycle_length = cycle_length
		self.neutron_source = neutron_source
		self.fuel_cycle = fuel_cycle 
		self.circuits = circuits
		self.target_burnup = burnup
		self.depl_parms = depl_parms
		self.pebble_flow = pebble_flow

		if MPO==None or cycle_length==None or neutron_source==None:
			poop

		self.mpos = self.mpo.split()

		# check to see that self.circuits and self.fuel_cycles match up

		total_volume, circ_volumes, prog_volumes = self.volume()

		for key in self.circuits:
			
			mpo1 = mocup.mpo()
			PF = {}
			for prog in self.fuel_cycle[key].keys():
				PF[prog] = self.pebble_flow[prog]
				for cell in self.fuel_cycle[key][prog]:
	 				mpo1 = self.mpos[cell] + mpo1

			circ = circuit()	
			circ.mixing = 'standard'
			circ.populate(fuel_cycle=self.fuel_cycle[key],MPO=mpo1, cycle_length=self.cycle_length*(circ_volumes[key]/total_volume), neutron_source=self.neutron_source, depl_parms=self.depl_parms, pebble_flow=PF)	
			self.passes[key] = circ

	def deplete(self):
		
		# This method depletes the passes in sequence

		import os
		
		for i in range(0,len(self.circuits)):

			# check to see if there is a BOP fuel vector

			ls = '\\rm tmp'
			print ls
			os.system(ls)

			ls = 'ls moi_files/. >> tmp'
			print ls
			os.system(ls)
			
			moi_files = open('tmp').read().split()

			if i == 0:
				for prog in self.fuel_cycle[self.circuits[0]].keys():
					ocf_loc = 'moi_files/moi.%s.0.pch' % self.fuel_cycle[self.circuits[0]][prog][0]
					bol = self.makeup*self.mpos[self.fuel_cycle[self.circuits[0]][prog][0]].volume[self.fuel_cycle[self.circuits[0]][prog][0]]
					bol.make_ocf(ocf_loc,lib=self.depl_parms.library)
	
			self.passes[self.circuits[i]].deplete()	

			if i + 1 < len(self.circuits):
				# There is an other pass that will need a BOC material
				#the discharge from the current circuit becomes the charge for the subsequence circuit

				# mix

				discharge = self.passes[self.circuits[i]].mix()

				# move materials correctly			

				for prog_key in self.fuel_cycle[self.circuits[i+1]].keys():
					BOP_cell = self.fuel_cycle[self.circuits[i+1]][prog_key][0]
					BOP_vol  = self.passes[self.circuits[i+1]].progressions[prog_key].depletions[BOP_cell].mpos[0].volume[BOP_cell]
					ocf_loc = 'moi_files/moi.%s.0.pch' % BOP_cell
					charge = discharge*BOP_vol
					charge.make_ocf(ocf_loc,lib=self.depl_parms.library)
			else:
				# there is a makeup fuel vector provided for the next circuit that will be created at the beginning of this loop or the discharge will not be reinserted into the core
				pass
	def power(self):
		
		# This method determines the power of an equilibrium fuel cycle by execiting the power method for each circuit in the fuel cycle
		
		self.ave_power = 0 

		for key in self.circuits:
			pow = self.passes[key].power()
			self.ave_power += pow
		return self.ave_power

	def burnup(self,p='no'):
		
		self.power()

		self.ave_BU = 0 
		total_volume = 0
		prog_volumes = {}
		circ_volumes = {}

		heavy_density = self.makeup.heavy_mass()*1e-6

		total_volume, circ_volumes, prog_volumes = self.volume()
	
		for circ_key in self.circuits:
			circ_burnup = 0
			circ_flow   = 0 
			prog_burnups = {}
			for prog_key in self.fuel_cycle[circ_key].keys():
				prog_burnups[prog_key] = 0
				for cell in self.fuel_cycle[circ_key][prog_key]:
					# note that since cycle lengths are now progression dependent I need to rederive this expression
					prog_burnups[prog_key] += self.passes[circ_key].progressions[prog_key].time[cell]*self.passes[circ_key].progressions[prog_key].ave_powers[cell]/(heavy_density*self.mpo.volume[cell])
				circ_burnup += prog_burnups[prog_key]*prog_volumes[circ_key][prog_key]*self.pebble_flow[prog_key]
				circ_flow   += prog_volumes[circ_key][prog_key]*self.pebble_flow[prog_key]
			self.ave_BU += circ_burnup/circ_flow

		if p == 'yes':
			output  = '\nCYCLE LENGTH SEARCH: PROGRESSION - %s, FLUX AMPLITUDE - %1.5E fission born neutrons per second, Target BURNUP %1.5E (MWd/MT)' % (self.title, self.neutron_source, self.target_burnup)
			output += '\nFUEL CYCLE INPUT PARMS:  CYCLE LENGTH - %1.5E EFPD, FLUX AMPLITUDE - %1.5E fission born neutrons per second' % (self.cycle_length, self.neutron_source)
			output += '\nFUEL CYCLE OUTPUT PARMS: POWER %1.5E MWth, BURNUP %1.5E MWd/MT' % (self.power(),self.ave_BU)
			open('scratch/DEPL.o','w').write(output)
			
		return self.ave_BU

	def normalized_burnup_estimate(self):

		# estimate power distribution using power density in mpos files (based on concentraitons from MCNP input deck and calculated cross sections)

		heavy_density = self.makeup.heavy_mass()*1e-6		   
		norm_BU = 0
		
		total_volume, circ_volumes, prog_volumes = self.volume()
		
		for circ_key in self.circuits:
			circ_burnup = 0
			prog_burnups = {}
			for prog_key in self.fuel_cycle[circ_key].keys():
				prog_burnups[prog_key] = 0
				for cell in self.fuel_cycle[circ_key][prog_key]:
					prog_burnups[prog_key] += self.mpos[cell].power_density[cell]*1.60217646e-19*self.neutron_source*self.mpos[cell].volume[cell]/(heavy_density*prog_volumes[circ_key][prog_key])*(circ_volumes[circ_key]/total_volume)
				circ_burnup += prog_burnups[prog_key]*prog_volumes[circ_key][prog_key]
			norm_BU += circ_burnup/circ_volumes[circ_key]
		return norm_BU

	def volume(self):
		
		total_volume = 0 
		circ_volumes = {}
		prog_volumes = {}
		
		for circ_key in self.circuits:
			circ_volumes[circ_key] = 0
			prog_volumes[circ_key] = {}
			for prog_key in self.fuel_cycle[circ_key].keys():
				prog_volumes[circ_key][prog_key] = 0
				for cell in self.fuel_cycle[circ_key][prog_key]:
					print self.mpos[cell].volume[cell]
					prog_volumes[circ_key][prog_key] += self.mpos[cell].volume[cell]
					circ_volumes[circ_key]           += self.mpos[cell].volume[cell]
					total_volume                     += self.mpos[cell].volume[cell]

		return total_volume, circ_volumes, prog_volumes

	def search_cycle_length(self,):

		# this module searches for the total cycle length that imposes a target burnup within error for a given target burnup,self.burnup , and source strength, self.neutorn_source

		import math
		import os

		os.system('\\rm scratch/cycle_length.o scratch/cycle_length.csv')

		# initialize cycle length output file
		os.system('touch scratch/cycle_length.csv')
		os.system('touch scratch/cycle_length.o')
		csv = '\ncycle length search\niteration (#),neutron source (n/s),power (MWt),cycle length (EFPD),burnup (MWd/MT)'
		open('scratch/cycle_length.csv','w').write(csv)		

		os.system('touch scratch/cycle_length.o')
		
		cycle_length = self.target_burnup/self.normalized_burnup_estimate()
	
		burnups = {(1.-self.depl_parms.burnup_perturb)*cycle_length:0,(1.+self.depl_parms.burnup_perturb)*cycle_length:0}
	
		k = 0
		for CL in sorted(burnups.keys()):
			k += 1
			self.populate(title=self.title, MPO=self.mpo, cycle_length=CL, neutron_source=self.neutron_source, fuel_cycle=self.fuel_cycle, circuits=self.circuits, depl_parms=self.depl_parms, burnup=self.target_burnup,pebble_flow=self.pebble_flow)
			self.deplete()
			burnups[CL] = self.burnup(p='yes')
			os.system('cat scratch/DEPL.o >> scratch/cycle_length.o')
			csv = '\n%d,%1.5E,%1.5E,%1.5E,%1.5E' % (k,self.neutron_source,self.power(),CL,burnups[CL])
			open('scratch/cycle_length.csv','a').write(csv)
 
		min = 0 
		max = 0 

		while math.fabs(self.burnup()/self.target_burnup-1) > self.depl_parms.burnup_err:

			k += 1
			
			m,b = least_squares(burnups)
			CL = (self.target_burnup - b)/m

			# ensure residence time is positive

			if CL < .01*sorted(burnups.keys())[0]:
				CL = .01*sorted(burnups.keys())[0]
			
			# ensure residence time is reasonably larger than the maximum cycle length estimate
			
			elif CL > 100*sorted(burnups.keys())[-1]:
				CL = 100.*sorted(burnups.keys())[-1]
			
			elif (min > 6) or (max > 6):
				CL = 0.5*(sorted(burnups.keys())[0] + sorted(burnups.keys())[-1])
			
			self.populate(title=self.title, MPO=self.mpo, cycle_length=CL, neutron_source=self.neutron_source, fuel_cycle=self.fuel_cycle, circuits=self.circuits, depl_parms=self.depl_parms, burnup=self.target_burnup,pebble_flow=self.pebble_flow)
			self.deplete()
			burnups[CL] = self.burnup(p='yes')

			os.system('cat scratch/DEPL.o >> scratch/cycle_length.o')
			csv = '\n%d,%1.5E,%1.5E,%1.5E,%1.5E' % (k,self.neutron_source,self.power(),CL,burnups[CL])
			open('scratch/cycle_length.csv','a').write(csv)

			if math.fabs(self.burnup()/self.target_burnup-1) < self.depl_parms.burnup_err:
				# cycle lenth has converged!
				pass
			else:
				# ensure best points for interpolation / extrapolation
				errs = {}
				for CL in burnups.keys():
					errs[math.fabs(burnups[CL] - self.target_burnup)] = CL	
				del burnups[errs[sorted(errs.keys())[-1]]]

				if self.burnup() > self.target_burnup:
					min  = 0 
					max += 1
				else:
					max  = 0 
					min += 1

		os.system('cat scratch/cycle_length.csv >> scratch/cycle_length.o')	

class circuit:

	# can be thought of as a pass through the core
	# a circuit is a collection of progressions 
	# each of these progressions has the same cycle length and the same charge fuel vector

	def __init__(self):
		self.progressions = {}
		self.cycle_length = 0 
		self.neutron_source = 0 

	def populate(self,fuel_cycle=None, MPO=None, cycle_length=None, neutron_source=None, depl_parms=None, pebble_flow=None):

		import mocup

		self.cycle_length = cycle_length
		self.neutron_source = neutron_source
		self.depl_parms = depl_parms

		if MPO==None or cycle_length==None or neutron_source==None:
			poop

		self.mpos = MPO.split()
		self.cycle_lengths = {}
		
		sum_flow = 0
		circ_vol = 0 	
		for prog_key in fuel_cycle.keys():
			prog_vol = 0 
			for cell in fuel_cycle[prog_key]:
				prog_vol += MPO.volume[cell]	
				circ_vol += MPO.volume[cell]
			sum_flow  += pebble_flow[prog_key]*prog_vol
		flow = sum_flow/circ_vol

		for key in fuel_cycle.keys():

			mpo1 = mocup.mpo()
			for cell in fuel_cycle[key]:
				mpo1 = self.mpos[cell] + mpo1

			self.cycle_lengths[key] = self.cycle_length*flow/pebble_flow[key]
			prog = progression()
			prog.populate(cells=fuel_cycle[key],MPO=mpo1, cycle_length=self.cycle_lengths[key], neutron_source=self.neutron_source, depl_parms=self.depl_parms)	
			self.progressions[key] = prog

	def deplete(self):

		
		for key in self.progressions.keys():

			self.progressions[key].deplete()

	def mix(self):
		
		import mocup

		total_volume = 0 
		mat = mocup.material()

		for key in self.progressions.keys():
			print key
			cell = self.progressions[key].cells[-1]
			total_volume += self.progressions[key].depletions[cell].mpos[0].volume[cell]

			mat1 = mocup.material()
			ocf_loc = 'moi_files/moi.%s.%d.pch' % (cell, (self.progressions[key].depletions[cell].mpos[-1].timestep+1))
			mat1.import_ocf(ocf_loc)

			mat = mat + mat1

		return mat*(1./total_volume)


	def power(self):
		
		# This method determines the power generate in this circuit by executing the power methods for each progression in the circuit

		self.ave_power = 0 
		
		for key in self.progressions.keys():
			pow = self.progressions[key].power()
			self.ave_power += pow

		return self.ave_power

class progression:
	def __init__(self):
		self.cells = []
		self.depletions = {}
	def populate(self, cells=None, MPO=None, cycle_length=None, neutron_source=None, depl_parms=None):

		import mocup
	
		if MPO==None or cells==None or neutron_source==None:
			poop

		self.cycle_length = cycle_length
		self.cells = cells
		self.depl_parms = depl_parms
		self.neutron_source = neutron_source
		self.mpos = MPO.split()
		self.time = {}
		total_volume = 0 

		for cell in self.cells:
			depl = mocup.Depletion()
			depl.mode    = self.depl_parms.mode
			depl.runTRANS = 'no'
			depl.method  = self.depl_parms.method
			depl.library = self.depl_parms.library
			depl.format  = self.depl_parms.format

			self.depletions[cell] = depl

			total_volume += self.mpos[cell].volume[cell]

		for cell in self.cells:

			self.time[cell] = self.cycle_length*self.mpos[cell].volume[cell]/total_volume
			
			# estimate the BOC mass in MT
			mat = mocup.material()
			mat.comp = self.mpos[cell].concentration[cell]
			mass = mat.heavy_mass()*self.mpos[cell].volume[cell]*(1e-6/.60221415)
			# estimate the expected power in MWt
			power = self.mpos[cell].power_density[cell]*self.mpos[cell].volume[cell]*self.neutron_source*1.6022e-19
			# estimate the burnup step in this cell
			expected_BU = (power*self.time[cell]/mass)
			
			if expected_BU < self.depl_parms.nBU1*self.depl_parms.dBU1 or expected_BU < self.depl_parms.dBU1*self.depl_parms.nBU:
				# expected burnup is less than the initial burnup steps
				# divide cycle length into self.nBU steps
				for i in range(self.depl_parms.nBU):
					self.depletions[cell].time.append(self.time[cell]/self.depl_parms.nBU)	
			else:
				# expected burnup is more than initial burnup steps (normal path)
				# generate timesteps based on not exceeding maximum burnup steps
				
				# define timesteps for intial burnup steps
				for i in range(self.depl_parms.nBU1):
					self.depletions[cell].time.append(self.depl_parms.dBU1*mass/power)

				expected_BU = expected_BU - self.depl_parms.dBU1*self.depl_parms.nBU1 
				nBU2 = int((expected_BU - expected_BU%self.depl_parms.dBU2)/self.depl_parms.dBU2 + 1.)
				if nBU2 < (self.depl_parms.nBU - self.depl_parms.nBU1):	
					nBU2 = int(self.depl_parms.nBU - self.depl_parms.nBU1)
				dBU2 = expected_BU/nBU2
				
				for i in range(nBU2):
					self.depletions[cell].time.append(dBU2*mass/power)
				
			for i in range(len(self.depletions[cell].time)):
				self.depletions[cell].source.append(self.neutron_source)
				self.depletions[cell].mpos.append(self.mpos[cell].twin())
				self.depletions[cell].mpos[i].timestep = i
				self.depletions[cell].start_ocf = self.depl_parms.start_ocf

	def deplete(self):

		import os
		import mocup

		#gofiss(self.depletions[self.cells[0]])

		origenBEAU(self.depletions[self.cells[0]])
		print 'gofiss %s' % self.cells[0]

		for i in range(1,len(self.cells)):
			
			mat = mocup.material()
			ocf_loc = 'moi_files/moi.%s.%d.pch' % (self.cells[i-1],self.depletions[self.cells[i-1]].mpos[-1].timestep+1)
			mat.import_ocf(ocf_loc)	
			
			mat = mat*(self.depletions[self.cells[i]].mpos[0].volume[self.cells[i]]/self.depletions[self.cells[i-1]].mpos[0].volume[self.cells[i-1]])

			ocf_loc = 'moi_files/moi.%s.0.pch' % (self.cells[i])
			mat.make_ocf(ocf_loc,lib=self.depl_parms.library)

			#gofiss(self.depletions[self.cells[i]])
			origenBEAU(self.depletions[self.cells[i]])

	def power(self):
	
		# This method analyzes the depletion objects and gives the time-averaged powers for each cell as well as the total time-averaged power for the progress
	
		import mocup
	
		# check to see that depletion has been performed

		for depl in self.depletions.values():
			if len(depl.mpos[0].ORIGEN_power.keys()) == 0:
				print 'deplete '
		
		self.ave_powers = {} 
		self.prog_power = 0 

		for cell in self.cells:

			self.ave_powers[cell] = 0 
			total_time = 0 

			if self.depl_parms.fission_power_only == 'yes':
				# use fission power only
				for i in range(len(self.depletions[cell].time)):
					total_time += self.depletions[cell].time[i]
					self.ave_powers[cell] += (self.depletions[cell].mpos[i].ORIGEN_fissionpower[cell]*self.depletions[cell].time[i])
			else:
				# use fission power and decay heat power
				for i in range(len(self.depletions[cell].time)):
					total_time += self.depletions[cell].time[i]
					self.ave_powers[cell] += (self.depletions[cell].mpos[i].ORIGEN_power[cell]*self.depletions[cell].time[i]) 

			self.ave_powers[cell]  = self.ave_powers[cell]/total_time
			self.prog_power       += self.ave_powers[cell]

		return self.prog_power

def gofiss(depletion):
	
	# deplete fuel in one cell of a fuel progression
	
	import mocup	
	
	for i in range(len(depletion.time)):
	
		# update the concentration vector in the mpos
		depletion.mpos[i].update()
		# run origen for sub state, i, in the depletion object
		mocup.origenBRO(depletion, i)
		

def origenBEAU(depletion):
	
	# This method depletes an entire irradiation history in one fell swoop
	#	This method should used in the case of constant flux / power AND constant cross sections or DECAY
	# 	This method reads the BOL material vector from a pch file, moi_files/moi.<cell>.0.pch

	import mocup	
	import math
	import os
	import re	

	nct = depletion.nct

	# import library directory
	dir = depletion.dir
	
	# import corss section libraries
	lib_loc = dir + depletion.library + '.lib'
	NXS_lib = open(lib_loc).read()
	k = NXS_lib.find('ACTINIDE')
	NXS_lib = NXS_lib[k:]
	
	# Link TAPE files to be used by origen

	# TAPE3 - blanket
	open('TAPE3.INP', 'w').write('')

	# TAPE9 - Neutron Cross Sections
	TAPE9 = 'cat %sdecay.lib %s%s.lib > TAPE9.INP' % (dir, dir, depletion.library)
	print TAPE9
	os.system(TAPE9)
	
	# TAPE10 - Gamma Cross Sectinos
	TAPE10 = 'ln -s %sgxuo2brm.lib TAPE10.INP' % (dir)
	print TAPE10
	os.system(TAPE10)
	
	mat = mocup.material()
	mat.lib = depletion.library
	mat.ORIGEN_nucl_list()
	
	for cell in depletion.mpos[0].cell:
	
		# intialize moi_file directory by deleting .out, .inp and .pch files for this cell

		tmp = '\\rm tmp'
		print tmp
		os.system(tmp)
		
		tmp = 'ls moi_files/moi.%s* >> tmp' % cell
		print tmp
		os.system(tmp)
		
		pch_loc = 'moi_files/moi.%s.0.pch' % cell
		for file in open('tmp').read().split():
			if file != pch_loc:
				remove = '\\rm %s' % file
				print remove
				os.system(remove)
	
		mpo = depletion.mpos[0]
		
		# Generate a list of nuclides in the same order as the ORIGEN cross section library
		
		mat = mocup.material()
		mat.lib = depletion.library
		mat.ORIGEN_nucl_list(lib=depletion.library)
		
		nuclides2 = []
		for nucl in depletion.mpos[0].mcnp_nuclides.values():
			nuclides2.append(depletion.nct[nucl])
		
		nuclides = []
	
		for nucl in (mat.structonly_nucl + mat.actin_nucl + mat.fisprod_nucl):
			if nucl in nuclides2:
				nuclides.append(nucl)
				if nucl in mpo.concentration[cell].keys():
					# There is a molar mass in the concentration vector, do nothing
					pass
				else:
					mpo.concentration[cell][depletion.nct[nucl]] = 0
		# identify which isotopes are not depleted so these can be added back into the material when the material vectors are updated
		mpo.structural_nuclides = []
		for nucl in depletion.mpos[0].mcnp_nuclides.values():
			if nucl not in nuclides:
				mpo.structural_nuclides.append(nucl)
		
		# update skeleton
		# initialize text for cross sections
		_SCROSS_SECTIONS_ = ''
		_ACROSS_SECTIONS_ = ''
		_FCROSS_SECTIONS_ = ''
		_STRUCTURAL_ = ''
		_ACTINIDES_ = ''
		_FISSIONPRODUCTS_ = ''
		
		# generate list of actinides and fission products to put in origen input deck
		
		for nucl in nuclides:
			moles = mpo.concentration[cell][nucl]*mpo.volume[cell]/.602214078
			if   nucl in mat.actin_nucl:
				_ACTINIDES_ += '        %s\n' % nucl
			elif nucl in mat.fisprod_nucl:
				_FISSIONPRODUCTS_ += '        %s\n' % nucl
			elif nucl in mat.structonly_nucl:
				_STRUCTURAL_ += '        %s\n' % nucl

		for nucl in nuclides:

			#determine if fuel or fission product and which library to append to and cross section format
			
			k = NXS_lib.split().index(nucl)
			if nucl in mat.actin_nucl:
				lib_num = mat.actinlib_num
			elif nucl in mat.fisprod_nucl:
				lib_num = mat.fisprodlib_num
			elif nucl in mat.structonly_nucl:
				lib_num = mat. structlib_num
			
			#determine metastable production ratio
			token = '      '
			token = token[:-len(nucl)] + nucl
			k = NXS_lib.find(token)
			
			# initialize fraction of NG XS and N2N XS that end in ground state and metastable states
			fNG    = 1
			fNGex  = 0
			fN2N   = 1
			fN2Nex = 0
			
			# intialize cross sections
			NG = 0 
			N2N = 0 
			NGex = 0 
			N2Nex = 0
			
			# grab cross sections from XS library to allocate the NG  and N2N cross section to reactions that end in ground state and metastable states
			if NXS_lib[k+14:k+16] != '  ':
				NG = float(NXS_lib[k+7:k+12])*math.pow(10,int(NXS_lib[k+13:k+16]))
			if NXS_lib[k+24:k+26] != '  ':
				N2N = float(NXS_lib[k+17:k+22])*math.pow(10,int(NXS_lib[k+23:k+26]))
			if NXS_lib[k+54:k+56] != '  ':
				NGex = float(NXS_lib[k+47:k+52])*math.pow(10,int(NXS_lib[k+53:k+56]))
			if NXS_lib[k+64:k+66] != '  ':
				N2Nex = float(NXS_lib[k+57:k+62])*math.pow(10,int(NXS_lib[k+63:k+66]))
			
			# determime fraction of reaction that end in ground state and metastable state
			if NG + NGex != 0:
				fNG   = NG/(NG + NGex)
				fNGex = NGex/(NG + NGex)
			
			if N2N + N2Nex != 0:
				fN2N   = N2N/(N2N+N2Nex)
				fN2Nex = N2Nex/(N2N+N2Nex)
			
			# Initialize Cross Sections at 0
			SN3N = 0
			SN2N = 0
			SNF  = 0
			SNG  = 0
			
			# define cross sections
			if nucl in mpo.N3N[cell].keys():
				SN3N = mpo.N3N[cell][nucl]
			if nucl in mpo.N2N[cell].keys():
				SN2N = mpo.N2N[cell][nucl]
			if nucl in mpo.NF[cell].keys():
				SNF = mpo.NF[cell][nucl]
			if nucl in mpo.NG[cell].keys():
				SNG = mpo.NG[cell][nucl]
			
			# write cross sections for a given nuclide
			
			# determine if a nuclide is actinide or fission product
			if nucl in mat.actin_nucl:
				# actinide
				#write line append to skeleton
				_ACROSS_SECTIONS_ += '%s %s %1.3E %1.3E %1.3E %1.3E %1.3E %1.3E -1.000000\n' % (lib_num, nucl, (fNG*SNG),(fN2N*SN2N),SN3N, SNF, (fNGex*SNG), (fN2Nex*SN2N))
				actinide_num = lib_num
			elif nucl in mat.fisprod_nucl:
				# fission produce
				token = '%s  %s' % (lib_num, nucl)
				k = NXS_lib.find(token)
				l = NXS_lib[k:].find('\n')
				m = NXS_lib[k+l+1:].find('\n')
				_XS_ = NXS_lib[k:(k+l+m+2)]
				# replace (N,G)
				SNGs = '%1.3E' % (SNG*fNG)
				_XS_ = _XS_[:12]+SNGs+_XS_[21:]
				SNGs = '%1.3E' % (SN2N*fN2N)
				_XS_ = _XS_[:22]+SNGs+_XS_[31:]
				SNGs = '%1.3E' % (SNG*fNGex)
				_XS_ = _XS_[:52]+SNGs+_XS_[61:]
				SNGs = '%1.3E' % (SN2N*fN2Nex)
				_XS_ = _XS_[:62]+SNGs+_XS_[71:]
				_FCROSS_SECTIONS_ += _XS_
				fissionproduct_num = lib_num
			elif nucl in mat.structonly_nucl:
				token = '%s        ' % (lib_num)
				token = token[:-len(nucl)] + nucl
				k = NXS_lib.find(token)
				l = NXS_lib[k:].find('\n')
				_XS_ = NXS_lib[k:(k+l+1)]
				SNGs = '%1.3E' % (SNG*fNG)
				_XS_ = _XS_[:12]+SNGs+_XS_[21:]
				SNGs = '%1.3E' % (SN2N*fN2N)
				_XS_ = _XS_[:22]+SNGs+_XS_[31:]
				SNGs = '%1.3E' % (SNG*fNGex)
				_XS_ = _XS_[:52]+SNGs+_XS_[61:]
				SNGs = '%1.3E' % (SN2N*fN2Nex)
				_XS_ = _XS_[:62]+SNGs+_XS_[71:]
				_SCROSS_SECTIONS_ += _XS_

		actinide_num = mat.actinlib_num
		fissionproduct_num = mat.fisprodlib_num
		structural_num = mat.structlib_num
	
		# combine cross sections from actinides and fission products
		_CROSS_SECTIONS_ = _SCROSS_SECTIONS_ + _ACROSS_SECTIONS_+_FCROSS_SECTIONS_
			
		# make ORIGEN input file(s)
			
		commands = 12
		irradiations = 0 
	
		header  =   ' -1'
		header += '\n -1'
		header += '\n -1'
		header += '\n TIT    %s' % depletion.title
		header += '\n BAS    BEAU point depletion calculation'
		header += '\n LIP    1 1 0 '
		if len(_STRUCTURAL_) > 0:
			header += '\n LPU \n' + _STRUCTURAL_ + '        -1'
		header += '\n LPU \n'
		header += _ACTINIDES_
		header += '        -1'
		header += '\n LPU \n' + _FISSIONPRODUCTS_ + '        -1'
		if len(_STRUCTURAL_) > 0:
			header += '\n LIB    0   1 2 3    -%s -%s -%s   9 50 0  4   0' % (structural_num, actinide_num, fissionproduct_num)
		else:
			header += '\n LIB    0   0 2 3    0 -%s -%s   9 50 0  4   0' % (actinide_num, fissionproduct_num)
		header += '\n OPTL   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
		header += '\n OPTA   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
		header += '\n OPTF   8 8 7 8 8  8 8 8 7 8  8 8 8 8 8  8 8 8 8 8  8 8 8 8'
		header += '\n CUT    3 1.0E-24  28 1.0E-75   -1'
		header += '\n INP    1   -2  -1  -1  1  1'
	
		deplet  = '\n _MODE_      _TIME1_   _POWER_    1   2   4   2'
		deplet += '\n _MODE_      _TIME2_   _POWER_    2   3   4   0'
		deplet += '\n _MODE_      _TIME3_   _POWER_    3   4   4   0'
		deplet += '\n _MODE_      _TIME4_   _POWER_    4   5   4   0'
		deplet += '\n _MODE_      _TIME5_   _POWER_    5   6   4   0'
		deplet += '\n _MODE_      _TIME6_   _POWER_    6   7   4   0'
		deplet += '\n _MODE_      _TIME7_   _POWER_    7   8   4   0'
		deplet += '\n _MODE_      _TIME8_   _POWER_    8   9   4   0'
		deplet += '\n _MODE_      _TIME9_   _POWER_    9  10   4   0'
		deplet += '\n _MODE_      _TIME10_   _POWER_   10  11   4   0'
		deplet += '\n OUT 11 1 0 0'
		deplet += '\n PCH 11 11 11'
		deplet += '\n MOV 11 1 0 1'

		footer  = '\n STP    4'
		footer += '\n' + _CROSS_SECTIONS_ + '0'

		TAPE5s = {}
	
		TAPE5 = header[:]
	
		BOL_index = 0 

		for i in range(len(depletion.time)):
			if commands + 13 < 300 and irradiations + 10 < 150:
				# the input doesn't violate the command limits
				pass
			else:
				# the input exceeds one of the command limits so a new input must be created
				commands = 12
				irradiations = 0 
				TAPE5 += footer
				TAPE5s[BOL_index] = TAPE5
				BOL_index = i

				# initialize new TAPE5
				TAPE5  = header[:]
				
			commands += 13
			irradiations += 10

			TAPE5 += deplet
			if depletion.mode == 'power':
				TAPE5 = TAPE5.replace('_MODE_','IRP')
				power_string = '%1.5E' % (depletion.power[i]*mpo.power[cell])
			elif depletion.mode == 'source' or depletion.mode == 'flux':
				TAPE5 = TAPE5.replace('_MODE_','IRF')
			 	if depletion.mode == 'source':
					power_string = '%1.5E' % (depletion.source[i]*mpo.flux_tally[cell])
				else:
					power_string = '%1.5E' % (depletion.power[i]*mpo.flux_tally[cell]/(mpo.global_power_density*1.60217464e-19))
			TAPE5 = TAPE5.replace('_POWER_',power_string)
			time_sum = 0 
			for j in range(1,11):
				time_sum    = j*depletion.time[i]/10
				token       = '_TIME%d_' % j
				time_string = '%1.5E' % time_sum
				TAPE5 = TAPE5.replace(token,time_string)

		TAPE5 += footer
		TAPE5s[BOL_index] = TAPE5
	
		# Begin depletion 
			
		# intialize inp and out files

		touch = 'touch moi_files/moi.%s.inp' % cell
		print touch
		os.system(touch)
			
		touch = 'touch moi_files/moi.%s.out' % cell
		print touch
		os.system(touch)
		
		for index in sorted(TAPE5s.keys()):
			
			# assume depletione.start_ocf == 0
			copy = '\cp moi_files/moi.%s.%d.pch TAPE4.INP' % (cell,index)
			print copy
			os.system(copy)
	
			open('TAPE5.INP','w').write(TAPE5s[index])
	
			# execute origen
				
			if depletion.library[0:3] in ['amo','emo','fft']:
				# origen uses a fast reactor set of cross sections, thus origen will be executed with o2_fast
				origen = 'o2_fast'
			else:
				# origen uses a thermal reactor set of cross sections, thus origen will be executed with o2_therm
				origen = 'o2_therm'
				
			print origen
			os.system(origen)

			TAPE6  = open( 'TAPE6.OUT').read()
			TAPE7  = open( 'TAPE7.OUT').read()			
			TAPE12 = open('TAPE12.OUT').read()	

			# parce TAPE6 for power and flux histories	
			ORIGEN_powers = re.compile('SP POW,MW\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+').findall(TAPE6)
			ORIGEN_fluxes = re.compile('NEUT. FLUX\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+').findall(TAPE6)

			# find locations of thermal heat summary tables
			tokens = re.compile('\s\d+\s{13,13}THERMAL').findall(TAPE12)
			decay_heat = []
			for token in tokens:
				skeleton_token = 'PAGE     '
				token1 = skeleton_token[:-len(token.split()[0])] + token.split()[0]
				k = TAPE6.find(token1)
				l = TAPE6[k:].find('0TOTAL')
				decay_heat.append( TAPE6[k+l:].split()[2:12])			

			# split TAPE7 into multiple PCH files, recover part that I 'split' out
			zeros = re.compile('\s{3,3}0\s{7,7}\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{5,5}[eE+-]{2,2}\d{2,2}\s+\d[.]\d{2,2}[eE+-]{2,2}\d{2,2}\s+').findall(TAPE7)
			pchs = TAPE7.split('0       0.00000E+00   0.00000E+00   0.00E+00')		

			i = 0 
			# organize data	
			for jndex in range(len(ORIGEN_powers)):

				if len(_STRUCTURAL_) > 0: 
					ORIGEN_power    = 0 
					ORIGEN_decayheat = 0 
					for j in range(3,len(ORIGEN_powers[jndex].split())): 
						ORIGEN_power += float(ORIGEN_powers[jndex].split()[j])/10.
						ORIGEN_decayheat += 1e-6*float(decay_heat[i][j-3])/10.
						ORIGEN_decayheat += 1e-6*float(decay_heat[i+1][j-3])/10.
						ORIGEN_decayheat += 1e-6*float(decay_heat[i+2][j-3])/10.

				else:
					ORIGEN_power    = 0
					ORIGEN_decayheat = 0
					for j in range(3,len(ORIGEN_powers[jndex].split())):
						ORIGEN_power += float(ORIGEN_powers[jndex].split()[j])/10.
						ORIGEN_decayheat += 1e-6*float(decay_heat[i  ][j-3])/10.
						ORIGEN_decayheat += 1e-6*float(decay_heat[i+1][j-3])/10.				
					i += 2

				ORIGEN_flux = 0
				for j in range(3,len(ORIGEN_fluxes[jndex].split())):
					ORIGEN_flux += float(ORIGEN_fluxes[jndex].split()[j])

				depletion.mpos[index+jndex].ORIGEN_power[cell]        = ORIGEN_power + ORIGEN_decayheat
				depletion.mpos[index+jndex].ORIGEN_fissionpower[cell] = ORIGEN_power
				depletion.mpos[index+jndex].ORIGEN_fluxes[cell]	      = ORIGEN_flux

				print (ORIGEN_power + ORIGEN_decayheat)

				# write EOS fuel composition				
				PCH = pchs[jndex] + zeros[jndex]
				pch_loc = 'moi_files/moi.%s.%d.pch' % (cell,(index+jndex+1))
				touch = 'touch %s' % pch_loc
				print touch
				os.system(touch)
				open(pch_loc,'w').write(PCH)
			
			# write ORIGEN output file to moi files
			move = 'cat TAPE12.OUT TAPE6.OUT >> moi_files/moi.%s.out' % cell
			print move
			os.system(move)

			# write ORIGEN input file to moi files			
			move = 'cat TAPE5.INP >> moi_files/moi.%s.inp' % cell	
			print move
			os.system(move)		

	# clean workspace
	remove = '\\rm TAPE*'
	print remove
	os.system(remove)

class AUSCE:
	
	def __init__(self):
		
		self.li6_tal   = 64 
		self.be9_tal   = 74
		self.li_mats   = []
		self.li6_ratio = []
		self.converged = 'no'
		self.mat_weighting = 0.50
		self.mat_err       = 0.01

	def equilibriate(self,depl_parms):
		
		# check convergence
		self.li6_convergence()

		if self.converged == 'no':
			
			import mocup
	
			if depl_parms.transport_module[0] in 'mM':
				output = open('outp1').read()
				input  = open('inp1').read()
				li6_tally, errs = self.fm4(output,self.li6_tal)
				be9_tally, errs = self.fm4(output,self.be9_tal)
				self.li6_ratio.append(be9_tally[be9_tally.keys()[0]][be9_tally[be9_tally.keys()[0]].keys()[0]]/li6_tally[li6_tally.keys()[0]][li6_tally[li6_tally.keys()[0]].keys()[0]])
				
			elif depl_parms.transport_module[0] in 'sS':
				pass
				# read input and output SERPENT deck
				# read reaction rate tallies

				
			for li_mat in self.li_mats:
				
				mat = mocup.material()
				if depl_parms.transport_module[0] in 'mM':
					mat.import_mcf('inp1',li_mat)
				elif depl_parms.transport_module[0] in 'sS':
					mat.import_scf('inp1',li_mat)

				if '40090' in mat.comp.keys() and '30060' in mat.comp.keys():
					token_mat = mat.mat_card

					if len(self.li6_ratio) > 1:
						li6_ratio = (1. - self.mat_weighting)*self.li6_ratio[-1] + self.mat_weighting*self.li6_ratio[-2]
					else:
						li6_ratio = self.li6_ratio[-1]
	
					mat.comp['30060'] = li6_ratio*mat.comp['40090']

					if depl_parms.transport_module[0] in 'mM':
						input = input.replace(token_mat,mat.mcf(li_mat))
					elif depl_parms.transport_module[0] in 'sS':
						input = input.replace(token_mat,mat.scf(li_mat))
	
			open('inp1','w').write(input)				

	def li6_convergence(self):
		
		from math import fabs
		if len(self.li6_ratio) >= 4:
			# enough iteratations
			if fabs(self.li6_ratio[-1]/self.li6_ratio[-2]-1) < self.mat_err:
				# li-6 ratio has converged
				self.converged = 'yes'

	def fm4(self,output,m):
		
		token = '1tally    '
		token = token[:-len(str(m))] + str(m)
		
		k = output.find(token)
		l = output[k:].find('\n\n multiplier bins') + len('\n\n multiplier bins')
		m = output[k+l:].find('\n\n')
		
		table_30 = output[k:k+l+m]
			
		cells = []
		bins = []
		
		k = table_30.find('cells') + len('cells')
		l = table_30[k:].find('\n\n')
		
		cells_string = table_30[k:k+l]
		
		while '(' in cells_string:
			k = cells_string.find('(')
			l = cells_string[k:].find(')')
			cell = cells_string[k:k+l+1]
			cell = cell.replace('\n','')
			while '  ' in cell:
				cell = cell.replace('  ',' ')
			cells.append(cell)
			cells_string = cells_string[:k] + cells_string[k+l+1:]
		
		cells += cells_string.split()

		if '  warning.' in table_30:
			k = table_30.find('  warning.')
			table_30 = table_30[:k]
		
		k = table_30.find('multiplier bins')

		for bin in table_30[k:].split('\n')[2:]:
			while '  ' in bin:
				bin = bin.replace('  ',' ')
			bins.append(bin[1:])
	
		while '' in bins:
			bins.remove('')

		k = output.find('1prob')
		l = output[k:].find(token)
		m = output[k+l:].find('==')
		
		tals = output[k+l:k+l+m].split('\n cell ')

		i = 1
		
		vals = {}
		errs = {}
		
		for cell in cells:
			vals[cell] = {}
			errs[cell] = {}
			for bin in bins:
				vals[cell][bin] = float(tals[i].split()[-2])
				errs[cell][bin] = float(tals[i].split()[-1])
				i += 1

		return vals, errs	
			

def search_comp(eq_fuel_cycle):

	# this method searches for the equilibrium state for a given equilbrium fuel cycle

	from os import system
	import mocup
	
	system('\\rm scratch/comp.o scratch/comp.csv')
	system('touch scratch/comp.o')
	system('touch scratch/comp.csv')
	csv = '\nequilibrium material search,\niterations (#),keff,convergence status'
	open('scratch/comp.csv','w').write(csv)

	if eq_fuel_cycle.initial_input == '':
		eq_fuel_cycle.initial_input == open('inp1').read()

	keff_list = []
	mat_list = []
	# this is a list to hold the trajectory of the equalibrium materials
	# mat_dict[<iteration>][<cell>] = mat

	# reset iteration counter and AUSCE parameters	
	i = 0 
	eq_fuel_cycle.AUSCE.converged = 'no'
	eq_fuel_cycle.AUSCE.li6_ratio = []	
	while 'Material Converged' not in comp_convergence(eq_fuel_cycle.depl_parms,keff_list,mat_list):
		i += 1
		csv = '\n%d,' % i
		open('scratch/comp.csv','a').write(csv)
		transportBRO(eq_fuel_cycle.depl_parms)
		mpo = mocup.mpo()
		mpo.populate(timestep=1)
		eq_fuel_cycle.mpo = mpo
		eq_fuel_cycle.search_source()
		eq_fuel_cycle.compBRO()
		eq_fuel_cycle.AUSCE.equilibriate(eq_fuel_cycle.depl_parms)
		keff_list.append(keff(eq_fuel_cycle.depl_parms))
		mat_list.append(eq_fuel_cycle.new_mats.copy())
		csv = '%1.5E,%s' % (keff(eq_fuel_cycle.depl_parms),comp_convergence(eq_fuel_cycle.depl_parms,keff_list,mat_list))
		open('scratch/comp.csv','a').write(csv)
		system('cat scratch/source.o >> scratch/comp.o')
	
	system('cat scratch/comp.csv >> scratch/comp.o')
	return keff(eq_fuel_cycle.depl_parms)

def comp_convergence(depl_parms, keff_list, mat_list):

	# this method checks the convergence of the the equilibrium cycle based on the number of iterations, convergences of keff and convergences of important nuclides in important cells.

	print depl_parms
	print dir(depl_parms)
	print depl_parms.import_nuclides.keys()

	from math import fabs
	
	if len(keff_list) > depl_parms.nMat and len(mat_list) > depl_parms.nMat:
		# enough iterations 
		if fabs((keff_list[-2] - keff_list[-1])/keff_list[-1]) < depl_parms.keff_err:
			# keffective has converged!!
	
			# check material composition
			for cell in depl_parms.import_nuclides.keys():
				for nucl in depl_parms.import_nuclides[cell]:
					if fabs((mat_list[-2][cell].comp[nucl] - mat_list[-1][cell].comp[nucl])/mat_list[-1][cell].comp[nucl]) > depl_parms.mat_err:
						return '%s in %s did not converge - %d: %1.5E, %d: %1.5E' % (nucl, cell, (len(mat_list)-2),mat_list[-2][cell].comp[nucl], (len(mat_list)-1), mat_list[-1][cell].comp[nucl])
			return 'Material Converged in %d iterations!'% len(keff_list)
		else:
			return 'keff has not converged - (n-1): %1.5E, (n): %1.5E' % (keff_list[-2],keff_list[-1])
	else:
		return 'not enough iterations - %d of %d iterations' % (len(keff_list), depl_parms.nMat)

def search_burnup(eq_fuel, burnup=None):
	
	# this method perturbs the burnup of the driver fuel 

	from math import fabs
	from os import system
	import mocup	

	# intialize scratch directory
	system('mkdir scratch')
	system('mkdir backup')
	system('mkdir moi_files')

	# initialize output files
	system('\\rm scratch/burnup.o')
	system('touch scratch/burnup.o')

	csv = 'burnup search\niterations (#),driver burnup (MWd/MT),eq. keff'
	open('scratch/burnup.csv','w').write(csv)

	# run transport to get mpo

	# import initital input deck
	if eq_fuel.initial_input == '':
		eq_fuel.initial_input = open('inp1').read()
	
	transportBRO(eq_fuel.depl_parms)
	mpo = mocup.mpo()
	if eq_fuel.depl_parms.transport_module[0] in 'mM':
		mpo.populate(timestep=1)
	elif eq_fuel.depl_parms.transport_module[0] in 'sS':
		mpo.populate_serpent(timestep=1)
	else:
		print 'undefined transport module!!'

	# search for the burnup of the driver circuit progress to impose a target power
	if burnup==None:
		# base target burnup of of the enrichment of the driver's makeup fuel
		driver_mat = eq_fuel.makeup[eq_fuel.driver]
		fissile = 0
		for nucl in ['922330','922350','942390','942410']:
			if nucl in driver_mat.comp.keys():
				fissile += driver_mat.comp[nucl]
		BU = fissile/driver_mat.heavy_moles()*1.e6
	else:
		BU = burnup
	keffs = {(1.-eq_fuel.depl_parms.burnup_perturb)*BU:0,(1.+eq_fuel.depl_parms.burnup_perturb)*BU:0,}

	# initialize a cycle length
	CL = {}
	for circ_prog_key in eq_fuel.fuel_cycle.keys():
		CL[circ_prog_key] = 1. # this is just a dummy value it will be over ridden by a reasonable estimate in the search_cycle_length loop 

	# initialize a neutron source 
	NS = eq_fuel.target_power/(200.*1.6022e-19)*2.5

	l = 0 
	for BU in keffs.keys():
		l += 1
		eq_fuel.target_burnups[eq_fuel.driver] = BU
		eq_fuel.populate(MPO=mpo, cycle_length=CL, neutron_source=NS)
		csv = '\n%d,%1.5E,' % (l,BU)
		open('scratch/burnup.csv','a').write(csv)
		keffs[BU] = search_comp(eq_fuel)
		system('cat scratch/comp.o >> scratch/burnup.o')
		csv = '%1.5E' % keffs[BU]
		open('scratch/burnup.csv','a').write(csv)

	min = 0 
	max = 0 

	while fabs(keff(eq_fuel.depl_parms) - eq_fuel.target_k) > eq_fuel.depl_parms.keff_err or l < 3.:
		
		# burnup has not converged
	
		l += 1
		m,b = least_squares(keffs)
		
		BU = (eq_fuel.target_k - b)/m

		# ensure burnup is positive
		
		if BU < .1*sorted(keffs.keys())[0]:	
			BU = .1*sorted(keffs.keys())[0]
		elif BU > 10*sorted(keffs.keys())[-1]:
			BU = 10*sorted(keffs.keys())[-1]
		elif (min > 6) or (max > 6):
			BU = 0.5*(sorted(keffs.keys())[0] + sorted(keffs.keys())[-1])
		
		eq_fuel.target_burnups[eq_fuel.driver] = BU
		eq_fuel.populate(MPO=mpo, cycle_length=CL, neutron_source=NS)	
		csv = '\n%d,%1.5E,' % (l,BU)
		open('scratch/burnup.csv','a').write(csv)
		keffs[eq_fuel.target_burnups[eq_fuel.driver]] = search_comp(eq_fuel)
		csv = '%1.5E' % keff(eq_fuel.depl_parms)
		open('scratch/burnup.csv','a').write(csv)
		system('cat scratch/comp.o >> scratch/burnup.o')
		
		if fabs(keff(eq_fuel.depl_parms) - eq_fuel.target_k) < eq_fuel.depl_parms.keff_err:
			
			# burnup has converged!!
			pass
	
		else: 
			# ensure best points for interpolation / extrapolation
			
			errs = {}
			for BU in keffs.keys():
				errs[fabs(keffs[BU] - eq_fuel.target_k)] = BU
			del keffs[errs[sorted(errs.keys())[-1]]]

			if keff(eq_fuel.depl_parms) > eq_fuel.target_k:
				min  = 0 
				max += 1
			else:
				max  = 0 
				min += 1
		
	output(eq_fuel)

def transportBRO(depl_parms, index=None):

	import os

	# this module runs neutron transport but also deletes current outp file

	if index == None:
		
		tmp = '\\rm tmp'
		print tmp
		os.system(tmp)
		
		tmp = 'touch tmp'
		print tmp
		os.system(tmp)
		
		tmp = 'ls backup/inp* >> tmp'
		print tmp
		os.system(tmp)
		
		scratch_indexes = open('tmp').read().replace(' ','').replace('\n','').split('backup/inp')
		if '' in scratch_indexes:
			i = scratch_indexes.index('')
			del scratch_indexes[i]
		if len(scratch_indexes) == 0:
			index = 0 
		else:
			for i in range(len(scratch_indexes)):
				scratch_indexes[i] = int(scratch_indexes[i])
			index = sorted(scratch_indexes)[-1]

	if depl_parms.runTRANS == 'yes':

		if depl_parms.transport_module[0] == 'm' or depl_parms.transport_module[0] == 'M':
			mcnp = '\\rm outp1 mctal1 runtp1' 
			print mcnp
			os.system(mcnp)
			mcnp = '\cp source srctp'
			print mcnp
			os.system(mcnp)
			mcnp = 'srun mcnp5.mpi i=inp1 o=outp1 mc=mctal1 runtpe=runtp1 >> mcnplog'
			print mcnp
			os.system(mcnp)
			mcnp = '\cp srctp source'
			print mcnp
			os.system(mcnp)
			mcnp = '\\rm runtp1 srctp' 
			print mcnp
			os.system(mcnp)

			# immediately move the output file into the backup directory
			mcnp = '\cp inp1 backup/inp%d' % (index + 1)
			print mcnp
			os.system(mcnp)
			
			mcnp = '\cp outp1 backup/outp%d' % (index + 1)
			print mcnp
			os.system(mcnp)		
		
		elif depl_parms.transport_module[0] == 's' or depl_parms.transport_module[0] == 'S':
			serpent = '\\rm inp1_res.m inp1.burn'
			print serpent
			#os.system(serpent)
			serpent = 'mpiexec sss117 inp1' 
			print serpent
			#os.system(serpent)
	
			# immediately move the output file into the backup directory
			serpent = '\cp inp1 backup/inp%d' % (index + 1)
			print serpent
			#os.system(serpent)
			
			serpent = '\cp inp1_res.m backup/inp%d_res.m' % (index + 1)
			print serpent
			#os.system(serpent)
			
			serpent = '\cm inp1.burn backup/inp%d.burn' % (index + 1)
			print serpent
			#os.system(serpent)

		else:
			print 'unknown transport module'
			exit()


def least_squares(y_vector):
        #generate sums

        sum_x = 0
        sum_y = 0
        sum_x2 = 0
        sum_xy = 0

        for x in y_vector.keys():

                sum_x += x
                sum_y += y_vector[x]
                sum_x2 += x**2
                sum_xy += x*y_vector[x]

        #generate slope and intercept

        slope = (sum_xy - 1./len(y_vector.values())*sum_x*sum_y)/(sum_x2 - 1./len(y_vector.values())*sum_x**2)
        x = y_vector.keys()[0]
        intercept = sum_y/len(y_vector.values()) - sum_x/len(y_vector.values())*slope

        return slope, intercept

def keff(depl_parm):

	if depl_parm.transport_module[0] in 'mM':	
		# this method returns the keff in the outp1 file
		outp = open('outp1',).read()
		k = outp.find('keff = ')
		return float(outp[k:].split()[2].replace(')',''))
	elif depl_parm.transport_module[0] in 'sS':
		outp = open('inp1_res.m')
		k = outp.find('IMP_KEFF')
		return float(outp[k:].split()[6])

def output(eq_fuel_cycle):
	
	# this file makes an output file of a given equilibrium fuel cycle

	import time

	output  = open('BEAU_art').read()
        output += '\n\n#initial time'
        output += "\ninitial_time = '"+eq_fuel_cycle.initial_time + "'"
        output += '\n\n#completion time'
        output += "\ncompletion_time = '" + time.ctime() + "'"


	output += '\n\n#'+'-'*79
	output += '\n# fuel cycle description'
	output += '\n#'+'-'*79

	output += '\n\n# target full core power (MWt)'
	output += '\ntarget_power = %1.5E' % eq_fuel_cycle.target_power

	output += '\n\n# target k'
	output += '\ntarget_k = %1.5E' % eq_fuel_cycle.target_k

	output += '\n\n# fuel cycle progression scheme'
	output +=   '\n#<circuit progression>'
	output +=   '\n#\t<circuit>'
	output +=   '\n#\t\t<fuel progression>'
	output +=   '\n#\t\t\t<burnup state cell ID>'

	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():	
		output += "\n{'" + circ_prog_key + "':\n\t{"
		for circ_key in eq_fuel_cycle.circuits[circ_prog_key]:
			output += "\n\t'" + circ_key + "':"
			output += '\n\t\t{'
			for prog_key in sorted(eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key].keys()):
				output += "\n\t\t'" + prog_key + "':\n\t\t\t["
				for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
					output += '%s,' % cell
				output += '],'
			output += '\n\t\t},'
		output += '\n\t},'
	output += '}'

	output += '\n\n#makeup material vectors'
	output += '\nimport mocup'
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		mat_ID = circ_prog_key.replace(' ','_')
		output += '\n%s = mocup.material()' % mat_ID 
		for nucl in eq_fuel_cycle.makeup[circ_prog_key].comp.keys():
			output += '\n%s.comp[%s] = %1.5E' % (mat_ID,nucl,eq_fuel_cycle.makeup[circ_prog_key].comp[nucl]) 

	output += '\n\n#'+'-'*79
	output += '\n# MCNP5 input deck'
	output += '\n#'+'-'*79
	output += '\n\n#\t' + eq_fuel_cycle.initial_input.replace('\n','\n#\t')

	output += '\n\n#'+'-'*79
	output += '\n# depletion parameters'
	output += '\n#'+'-'*79

	output += '\n\n# minimum number of depletion burnup steps'
	output += '\nnBU = %d' % eq_fuel_cycle.depl_parms.nBU

	output += '\n\n# minimum number of initial depletion burnup steps'
	output += '\nnBU1 = %d' % eq_fuel_cycle.depl_parms.nBU1

	output += '\n\n# initial burnup steps (MWd/MT)'
	output += '\ndBU1 = %1.5E' % eq_fuel_cycle.depl_parms.dBU1
	
	output += '\n\n# maximum burnup step (MWd/MT)'
	output += '\ndBU2 = %1.5E' % eq_fuel_cycle.depl_parms.dBU2
	
	output += '\n\n# auxillary cross section library'
	output += '\n#\tthese cross sections are used for all the reactions without system specific cross sections'
	output += "\nlibrary = '%s'" % eq_fuel_cycle.depl_parms.library
	
	output += '\n\n# origen flux mode'
	output += '\n#\tsource = constant flux based on estimate of the neutron fission source (recommended)'
	output += '\n#\tflux   = constant flux calculated based on power density from mpo to impose a target power level'
	output += '\n#\tpower  = constant power calculated based on power distribution from mpo'
	output += "\nmode = '%s'" % eq_fuel_cycle.depl_parms.mode

	output += '\n\n# origen depletion step methoding'
	output += '\n#\tbeginning = BOS material depleted based on cross section calculated at BOS (reccommened bc BUS cross sections are over ridden onto mpos)'
	output += '\n#\tpredictor corrector = MOS material is generated to generate MOS cross sections, these cross sections are used to deplete fuel from BOS to EOS'
	output += "\nmethod = '%s'" % eq_fuel_cycle.depl_parms.method
	

	output += '\n\n# flag to write or not write mpos during the Equilibrium Fuel Cycle Search'
	output += "\nfortmat = '%s'" % eq_fuel_cycle.depl_parms.format
	
	output += '\n\n# transport module'
	output += "\ntransport_modulue = '%s'" % eq_fuel_cycle.depl_parms.transport_module

	output += '\n\n#the error acceptable in terms of the target burnup'
	output += '\nburnup_err = %1.5E' % eq_fuel_cycle.depl_parms.burnup_err
	
	output += '\n\n# the amount by which BEAU perturbs cycle length and burnup in search_cyclelength() and search_burnup() modules'
	output += '\nburnup_perturb = %1.5E' % eq_fuel_cycle.depl_parms.burnup_perturb

	output += '\n\n# the error acceptable in terms of the target power'
	output += '\npower_err = %1.5E' % eq_fuel_cycle.depl_parms.power_err
	
	output += '\n\n# the amount by which BEAU perturbs the neutron source intensity in search_source() module'
	output += '\npower_perturb = %1.5E' % eq_fuel_cycle.depl_parms.power_perturb
	
	output += '\n\n# the error acceptable in isotopic concentration in terms of the best estimate of isotopic concentration for an important isotope in a checked cell'
	output += '\nmat_err = %1.5E' % eq_fuel_cycle.depl_parms.mat_err
	
	output += '\n\n# the error acceptable in k-effective'
	output += '\nkeff_eff = %1.5E' % eq_fuel_cycle.depl_parms.keff_err

	output += '\n\n# minimum number of guesses at the equilibrium state before BEAU is allowed to converge on an equilibrium state'
	output += '\nnMat = %d' % eq_fuel_cycle.depl_parms.nMat
	
	output += '\n\n# material weighting factor'
	output += '\n#\tthis is the fraction by which the i+1 equilibrium compositions are weighted by the i equilibrium composition'
	output += '\nmat_weighting = %1.5E' % eq_fuel_cycle.depl_parms.mat_weighting
	
	output += '\n\n# important nuclides dictionary'
	output += '\n#\tthis lists the nuclides in which cells will be tested for convergence'
	output += '\n#<cell>'
	output += '\n#\t<isotope>'
	output += '\nimport_nuclides = {'
	for cell in eq_fuel_cycle.depl_parms.import_nuclides.keys():
		output += "\n'" + cell + "':[\n\t"
		for nucl in eq_fuel_cycle.depl_parms.import_nuclides[cell]:
			output += "'" + nucl + "',"
		output += '],'
	output += '\n}'

	output += '\n\n#tally ID for li6 (n,abs) reaction'
	output += '\nli6_tal = %d' % (eq_fuel_cycle.AUSCE.li6_tal)
	
	output += '\n\n#tally ID for be9 (n,alpha) reaction'
	output += '\nbe9_tal = %d' % eq_fuel_cycle.AUSCE.be9_tal

	output += '\n\n# materials that BEAU sets to the calculated li6 ratio'
	output += '\nli_mats = ['
	for mat_ID in eq_fuel_cycle.AUSCE.li_mats:
		output += '%d,' % int(mat_ID)
	output += ']'
	
	output += '\n\n# lithium material weighting factor'
	output += '\n#\tthis is the fraction by which the i+1 li ratio is weighted by the i li ratio'
	output += 'li_weighting = %1.5E' % eq_fuel_cycle.AUSCE.mat_weighting
	
	output += '\n\n# the error acceptable in terms of the best estimate of the li ratio'
	output += '\nli_err = %1.5E' % eq_fuel_cycle.AUSCE.mat_err

	output += '\n\n#'+'-'*79
	output += '\n# fuel cycle search results'
	output += '\n#'+'-'*79

	output += '\n\n# Equilibrium k-effective'
	output += 'eq_keff = %1.5E' % keff(eq_fuel_cycle.depl_parms)
	
	output += '\n\n# full core power (MWt)'
	output += '\npower = %1.5E' % eq_fuel_cycle.power()	
	output += '\n\n# power distribution'
	output += '\n#<cell>\n#\t<power (MWt>\npower_distribution = {'
	
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		for circ_key in eq_fuel_cycle.circuits[circ_prog_key]:
			for prog_key in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key].keys():
				for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
					output += "\n'" + cell + "':\n\t"
					output += '%1.5E,' % eq_fuel_cycle.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].ave_powers[cell]
	output += '\n}'

	output += '\n\n# neutron fission source strength (n/s)'
	output += '\nsource = %1.5E' % eq_fuel_cycle.neutron_source
	output += '\n\n# flux distrubiton (n/cm2s)\n#<cell>\n#\t<flux>'
	output += '\nflux_distribuiton = {'
	
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		for circ_key in eq_fuel_cycle.circuits[circ_prog_key]:
			for prog_key in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key].keys():
				for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
					output += "\n'" + cell + "':\n\t"
					output += '%1.5E,' % (eq_fuel_cycle.neutron_source*eq_fuel_cycle.mpo.flux_tally[cell])
	output += '\n}'	

	output += '\n\n# burnups (MWd/MT)'
	output += '\n#<circuit progression>\n#\t<burnup>'
	output += '\nburnup = {'
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		output += "'%s':%1.5E," % (circ_prog_key,eq_fuel_cycle.circ_progs[circ_prog_key].burnup())
	output += '\n}'
	output += '\n\n#burnup distribution (fima %)'
	output += '\n#<cell>\n#\t<burnup>'
	output += '\nfima = {'
	
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		a = eq_fuel_cycle.makeup[circ_prog_key]*.60221415
		for circ_key in eq_fuel_cycle.circuits[circ_prog_key]:
			for prog_key in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key].keys():
				for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
					output += "\n'" + cell + "':\n\t"
					output += '%2.1f,' % (100*(1.-eq_fuel_cycle.new_mats[cell].heavy_moles()/a.heavy_moles()))
	output += '\n}'
	
	output += '\n\n# cycle lengths (EFPD)'
	output += '\n#<circuit progression>\n#\t<burnup>'
	output += '\ncycle_lengths = {'
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		output += "'%s':%1.5E" % (circ_prog_key,eq_fuel_cycle.cycle_length[circ_prog_key]) 
	output += '\n}'
	output += '\n\n#residence time in burnup states (EFPD)'
	output += '\n#<cell>\n#\t<burnup>'
	output += '\nresidence_times = {'
	
	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		for circ_key in eq_fuel_cycle.circuits[circ_prog_key]:
			for prog_key in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key].keys():
				for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
					output += "\n'" + cell + "':\n\t"
					output += '%1.5E,' % eq_fuel_cycle.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].time[cell]
	output += '\n}'

	open('output.py','w').write(output)

	import mocup

	mat = mocup.material()
	mat.ORIGEN_nucl_list()

	for circ_prog_key in eq_fuel_cycle.fuel_cycle.keys():
		for circ_key in eq_fuel_cycle.circuits[circ_prog_key]:
			for prog_key in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key].keys():
				mat_out  = "This is the material concentration vector output file for progression '%s'" % prog_key
				mat_out += '\n\ncells:,'
				heading1 = '\n ,'
				heading2 = '\nisotope,'
				sum_time = 0 
				for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
					heading1 += 'burnup state: %s,' % cell
					mat_out += '%s,' % cell
					heading2 += '%1.5E,' % sum_time
					for dt in eq_fuel_cycle.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].time:
						sum_time += dt
						heading1 += ' ,'
						heading2 += '%1.5E,' % sum_time

				mat_out += '\n' + heading1 + heading2

				nuclides = eq_fuel_cycle.mpo.concentration[cell].keys()
				mat.ORIGEN_nucl_list()
				for nucl in nuclides:
					if nucl not in eq_fuel_cycle.mpo.NG[cell].keys() or nucl not in (mat.actin_nucl + mat.fisprod_nucl):
						print nucl
						nuclides.remove(nucl)
	
				for nucl in nuclides:
					mat_out += '\n%s,' % nucl
					for cell in eq_fuel_cycle.fuel_cycle[circ_prog_key][circ_key][prog_key]:
						for i in range(len(eq_fuel_cycle.circ_progs[circ_prog_key].passes[circ_key].progressions[prog_key].depletions[cell].time)+1):
							pch_loc = 'moi_files/moi.%s.%d.pch' % (cell,i)
							mat = mocup.material()
							mat.import_ocf(pch_loc)
							mat = mat*(.60221415/eq_fuel_cycle.mpo.volume[cell])	
							if nucl in mat.comp.keys():
								mat_out += '%1.5E,' % mat.comp[nucl]
							else:
								mat_out += '%1.5E,' % 0

				mat_loc = 'concentration.%s.%s.%s.csv' % (circ_prog_key,circ_key,prog_key)

				open(mat_loc,'w').write(mat_out)




