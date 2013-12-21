#!/usr/bin/env python

import BEAU
import mocup
import os

# -----------------------------------------------------------------------------
# DEFINE EQUILIBRIUM DEPLETION ANALYSIS PROBLEM
# -----------------------------------------------------------------------------

# initialize Equilibrium Depletion
Eq_Cycle = BEAU.equilibrium_cycle()

# DEFINE FUEL PROGRESSION
#<circulation progression>
#	<circulation>
#		<progression>
#			<list of burnup state>
Eq_Cycle.fuel_cycle  = { 
'LEU':{
	'pass1':{
		'progression11':
			['11100','12100','13100','14100','15100'],
		'progression21':
			['21100','22100','23100','24100','25100'],
		'progression31':
			['31100','32100','33100','34100','35100'],
		'progression41':
			['41100','42100','43100','44100','45100'],
		},
	'pass2':{
		'progression12':
			['11200','12200','13200','14200','15200'],
		'progression22':
			['21200','22200','23200','24200','25200'],
		'progression32':
			['31200','32200','33200','34200','35200'],
		'progression42':
			['41200','42200','43200','44200','45200'],
		},
        'pass3':{
                'progression13':
                        ['11300','12300','13300','14300','15300'],
                'progression23':
                        ['21300','22300','23300','24300','25300'],
                'progression33':
                        ['31300','32300','33300','34300','35300'],
                'progression43':
                        ['41300','42300','43300','44300','45300'],
                },
        'pass4':{
                'progression14':
                        ['11400','12400','13400','14400','15400'],
                'progression24':
                        ['21400','22400','23400','24400','25400'],
                'progression34':
                        ['31400','32400','33400','34400','35400'],
                'progression44':
                        ['41400','42400','43400','44400','45400'],
                },
        'pass5':{
                'progression15':
                        ['11500','12500','13500','14500','15500'],
                'progression25':
                        ['21500','22500','23500','24500','25500'],
                'progression35':
                        ['31500','32500','33500','34500','35500'],
                'progression45':
                        ['41500','42500','43500','44500','45500'],
                },
        'pass6':{
                'progression16':
                        ['11600','12600','13600','14600','15600'],
                'progression26':
                        ['21600','22600','23600','24600','25600'],
                'progression36':
                        ['31600','32600','33600','34600','35600'],
                'progression46':
                        ['41600','42600','43600','44600','45600'],
                },
        'pass7':{
                'progression17':
                        ['11700','12700','13700','14700','15700'],
                'progression27':
                        ['21700','22700','23700','24700','25700'],
                'progression37':
                        ['31700','32700','33700','34700','35700'],
                'progression47':
                        ['41700','42700','43700','44700','45700'],
                },
        'pass8':{
                'progression18':
                        ['11800','12800','13800','14800','15800'],
                'progression28':
                        ['21800','22800','23800','24800','25800'],
                'progression38':
                        ['31800','32800','33800','34800','35800'],
                'progression48':
                        ['41800','42800','43800','44800','45800'],
                },
	}
}

Eq_Cycle.target_burnups = {'LEU':0,}
Eq_Cycle.target_power = 2.36000E+02

# define the order of circulations
Eq_Cycle.circuits = {'LEU':['pass1','pass2','pass3','pass4','pass5','pass6','pass7','pass8']}
	
# define the relative pebble flux (pebble/(unit area unit time))
#	if not defined a flat distribution will be assumed
Eq_Cycle.pebble_flow = {
        'progression11':1.,
        'progression21':1.,
        'progression31':1.,
        'progression41':1.,
        'progression12':1.,
        'progression22':1.,
        'progression32':1.,
        'progression42':1.,
        'progression13':1.,
        'progression23':1.,
        'progression33':1.,
        'progression43':1.,
        'progression14':1.,
        'progression24':1.,
        'progression34':1.,
        'progression44':1.,
	'progression15':1.,
	'progression25':1.,
	'progression35':1.,
	'progression45':1.,
	'progression16':1.,
	'progression26':1.,
	'progression36':1.,
	'progression46':1.,
        'progression17':1.,
        'progression27':1.,
        'progression37':1.,
        'progression47':1.,
        'progression18':1.,
        'progression28':1.,
        'progression38':1.,
        'progression48':1.,
}

# define the makeup fuel vector for each circulation progression

# initiate material object
seed = mocup.material()
# define its composition vector

seed.comp['922350'] = 7.89393E-03
seed.comp['922380'] = 3.13735E-02
seed.comp['541350'] = 1e-30


# link material object to the progression key
Eq_Cycle.makeup = {'LEU':seed,}

# define the ID of the driver circuit progression

Eq_Cycle.driver = 'LEU'

# DEFINE THE STANDARD ORIGEN LIBRARY TO BE USED
Eq_Cycle.depl_parms.library = 'pwru50'
Eq_Cycle.depl_parms.dBU1 = 1000
Eq_Cycle.depl_parms.dBU2 = 10000
Eq_Cycle.depl_parms.nBU1 = 3
Eq_Cycle.depl_parms.import_nuclides = {
        '13100':['922350','922380','942390'],
        '13200':['922350','922380','942390'],
        '13300':['922350','922380','942390'],
        '13400':['922350','922380','942390'],
	'13500':['922350','922380','942390'],
	'13600':['922350','922380','942390'],
	'13700':['922350','922380','942390'],
	'13800':['922350','922380','942390'],
}
Eq_Cycle.depl_parms.mat_err = 0.1

Eq_Cycle.AUSCE.li_mats = [5,6,23,24,25,700,]

BEAU.search_burnup(Eq_Cycle,burnup=1.77878E+05)



