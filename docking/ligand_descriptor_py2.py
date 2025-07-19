#! /usr/bin/env python2.4

###################################################################################

"""
Aggiungere funzione per alcolare stato di protonazione
"""

def PainsFinder( smile, smart_def_file ):

        import sys
        sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
        import openbabel, pybel

	pymol = pybel.readstring( "smi", smile )
	pymol.removeh()

	pains_match = 'No'
	substr_name = []
	highlight   = []
	matches     = []

 	input = open( smart_def_file,'r' )

	for line in input.readlines():
		data = line.split()
		smarts = pybel.Smarts( data[0] )
		matches = smarts.findall( pymol )

		if matches:
			pains_match = 'Yes'
			atoms       = []
			info        = ''

			for i in range( 0, len(matches) ):
				for j in range( 0, len(matches[i]) ):
					atoms.append( matches[i][j] )
				atoms.sort()
			highlight.append( atoms )

			for i in range( 1, len(data) ):
				if i == 1:
					info = data[i]
				else:
					info = info + ' ' + data[i]
			substr_name.append( info )

	input.close() 

        return pains_match, substr_name, highlight

###################################################################################

def GetInChI( smile ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )

	import openbabel

        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("smi", "inchi")

        mol = openbabel.OBMol()
        conv.ReadString(mol,smile)

        inchi = conv.WriteString(mol)

	inchi = inchi.rstrip( '\n' )

        return inchi

###################################################################################

def GetInChIKey( smile ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )

	import openbabel

        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("smi", "inchikey")

        mol = openbabel.OBMol()
        conv.ReadString(mol,smile)

        inchikey = conv.WriteString(mol)

	inchikey = inchikey.rstrip( '\n' )

        return inchikey

###################################################################################

def ChiralAtoms( OBMol ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel
	
	n = 0
	chiral = OBMol.IsChiral()	
	if chiral:
		for obatom in openbabel.OBMolAtomIter( OBMol ):
			if obatom.IsChiral() and obatom.GetType() != 'N3':
				n = n + 1

	return chiral, n

###################################################################################

def RingsCount( OBMol ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	n_rings           = 0
	n_aromatic_rings  = 0
	n_aliphatic_rings = 0

	for ring in openbabel.OBMolRingIter( obmol ):
		n_rings = n_rings + 1
		if ring.IsAromatic():
			n_aromatic_rings = n_aromatic_rings + 1 
		else:
			n_aliphatic_rings = n_aliphatic_rings + 1

	return n_rings, n_aromatic_rings, n_aliphatic_rings
	
###################################################################################

def logBB( logP, tpsa ):

	'''
	Vilar, S.; Chakrabarti, M.; Costanzi, S.
	Prediction of passive blood-brain partitioning: straightforward and 
	effective classification models based on in silico derived 
	physicochemical descriptors.
	J Mol Graph Model, 2010, 8, 899-903. DOI:10.1016/j.jmgm.2010.03.010

	CNS availability: model_1 >= 0.3 | model_2 >= 1
	'''

	model_1 = 0.0
	model_1 = 0.5159 * logP - 0.0277 * tpsa - 0.3462 

	'''
	model_2 = 0.0
	model_2 = 0.2289 * logP - 0.5671 * ( n_acid_groups + n_base_group ) + 2.3420
	'''

	return model_1

###################################################################################

def Lipinski( HBA, HBD, logP, MW ): 

	''' 

	Lipinski, C.A.; Lombardo, F.; Dominy, B.W.; Feeney, P.J. 
	Experimental and computational approaches to estimate solubility and 
	permeability in drug discovery and development settings. 
	Adv Drug Del Rev, 1997, 46, 3-26 - DOI:10.1016/S0169-409X(00)00129-0

	'''    

	n = 0
	rule_of_3 = False
	rule_of_5 = False

	if HBA <= 5:
		n = n + 1
	if HBD <= 10:
		n = n + 1
	if logP <= 5.0:
		n = n + 1
	if MW <= 500.0:
		n = n + 1

	if n == 4:
		rule_of_3 = True
		rule_of_5 = True
	if n == 3:
		rule_of_3 = True
		rule_of_5 = False

	return rule_of_3, rule_of_5

###################################################################################

def IsChiralitySpecified( smile, n_chiral_atms ):

	defined_atms = 0
	def_pos = ''

	if n_chiral_atms < 1:
		specified = False
		return specified

	for n in range( 0, len( smile ) ):
		cha = smile[n]
		if cha == '@':
			defined_atms = defined_atms + 1
			if def_pos == ( n - 1 ):
				defined_atms = defined_atms - 1
			def_pos = n

	if defined_atms == n_chiral_atms:
		specified = True
	else:
		specified = False

	return specified

###################################################################################

def DoubleBondsAna( OBMol, d_bonds, smile ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	d_bonds_in_ring = 0
	d_bonds_out_ring = 0
	defined_EZ = False
	n_def_signs = 0

	if d_bonds > 0:
		for bond in openbabel.OBMolBondIter( OBMol ):
			if bond.IsDouble():
				if bond.IsInRing():
					d_bonds_in_ring = d_bonds_in_ring + 1
				else:
					sa = ( bond.GetBeginAtom() ).IsCarbon()
					fa = ( bond.GetEndAtom() ).IsCarbon()
					if sa and fa:
						d_bonds_out_ring = d_bonds_out_ring + 1

	if d_bonds_out_ring > 0:
		for n in range( 0, len( smile ) ):
			cha = smile[n]
			if cha == '/' or cha == '\\':
				n_def_signs = n_def_signs + 1
		if n_def_signs == ( d_bonds_out_ring * 2 ) or n_def_signs - 1 == ( d_bonds_out_ring * 2 ) :
			defined_EZ = True

	return d_bonds_in_ring, d_bonds_out_ring, defined_EZ

	
###################################################################################

def my_IsChetonGroup( bond ):

        import sys 
        sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
        import openbabel

        atomA     = bond.GetBeginAtom()
        atomTypeA = atomA.GetType()

        atomB     = bond.GetEndAtom()
        atomTypeB = atomB.GetType()

        if atomTypeA == 'C2':
                carbon = atomA
        else:
                carbon = atomB

	if carbon.GetHvyValence() < 3:
		return False

	if carbon.GetHvyValence() == 3:
		for neighbour_atom in openbabel.OBAtomAtomIter( carbon ):
			if neighbour_atom.GetType() == 'O3':
                		return False
	
        return True

###################################################################################

def my_IsAldehydeGroup( bond ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	atomA     = bond.GetBeginAtom()
	atomTypeA = atomA.GetType()

	atomB     = bond.GetEndAtom()
	atomTypeB = atomB.GetType()

	if atomTypeA == 'C2':
		carbon = atomA
	else:
		carbon = atomB

	if carbon.GetHvyValence() <= 2:
		return True
	else:
		return False

###################################################################################

def my_IsNotCarboxylGroup( bond ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	atomA     = bond.GetBeginAtom()
	atomTypeA = atomA.GetType()

	atomB     = bond.GetEndAtom()
	atomTypeB = atomB.GetType()

	if atomTypeA == 'O2':
		oxygen = atomA
	else:
		oxygen = atomB

	if not oxygen.IsCarboxylOxygen():
		return True
	else:
		return False

###################################################################################

def my_GetFingerPrint( OBMol ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	obConversion.SetOutFormat( "fpt" )
	obConversion.SetOptions( "xs", obConversion.OUTOPTIONS )
	fpt = obConversion.WriteString( OBMol )

	return fpt[1:]

###################################################################################

def format_A_to_format_B( format_in, format_out, mol_in ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats( format_in, format_out )

	OBMol_out = openbabel.OBMol()
	obConversion.ReadString( OBMol_out, mol_in )

	mol_out = obConversion.WriteString( OBMol_out )

	return mol_out

###################################################################################

def my_3DBuilder( format, smile ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel
	
	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats( "smi", format )

	OBMol = openbabel.OBMol()
	obConversion.ReadString( OBMol, smile )

	bd = openbabel.OBBuilder()
	bd.Build( OBMol )

	OBMol.SetDimension( 3 )
	OBMol.AddHydrogens()

	ff = openbabel.OBForceField.FindForceField( "MMFF94" )
	ff.Setup( OBMol )
	ff.SteepestDescent( 150, 1.0e-4 )
	ff.WeightedRotorSearch( 5, 25 )
 	ff.ConjugateGradients( 250, 1.0e-6 )
	ff.UpdateCoordinates( OBMol )
	
	OBMol.SetTitle( "CMLDID" )
	mol3D = obConversion.WriteString( OBMol )

	return mol3D

###################################################################################

def my_2DBuilder( format, smile ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	obConversion = openbabel.OBConversion()
	obConversion.SetInAndOutFormats( "smi", format )

	OBMol = openbabel.OBMol()
	obConversion.ReadString( OBMol, smile )

	openbabel.OBOp.gen2d = openbabel.OBOp.FindType( "Gen2D" )
	openbabel.OBOp.gen2d.Do( OBMol )

	OBMol.DeleteHydrogens()
	OBMol.SetTitle( "CMLDID" )

	mol2D = obConversion.WriteString( OBMol )

	return mol2D

###################################################################################

def IsComplex( smile ):

	import sys
	sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )
	import openbabel

	obConversion = openbabel.OBConversion()
	obConversion.SetInFormat("smi")

	obmol = openbabel.OBMol()
	obConversion.ReadString( obmol, smile )

	n_ha_atms  = obmol.NumHvyAtoms()
	n_ha_rb    = obmol.NumRotors()

	if n_ha_atms > 15 or n_ha_rb > 5:
		complexity = True
	else:
		complexity = False

	return complexity

###################################################################################

def check_web_data( web_user_rightconfig_spec, web_user_delivery, chiral, specified_chirality, d_bonds_out_ring, defined_EZ ):

	import sys

	if web_user_rightconfig_spec == "NotChiral" and d_bonds_out_ring > 0.0:
		print "ERROR| Double bond(s) detected but \"No isomers\" has been declared"
		sys.exit()

	if web_user_rightconfig_spec == "NotChiral" and d_bonds_out_ring > 0.0 and defined_EZ == False:
		print "ERROR| Double bond(s) detected but geometry has been not specified in molecule structure"
		sys.exit()

	if web_user_rightconfig_spec == "NotChiral" and chiral == True:
		print "ERROR| Chiral atom(s) detected but \"No isomers\" has been declared"
		sys.exit()

	if web_user_delivery == "Theoretical":
		if chiral == True and specified_chirality == False:
			print "ERROR| Chirality must be specified for theoretical compounds"
			sys.exit()
		if d_bonds_out_ring > 0.0 and defined_EZ == False:
			print "ERROR| Double bonds geometry must be specified for theoretical compounds"
			sys.exit()
		if web_user_rightconfig_spec == "No":
			print "ERROR| Theoretical entry must be declared (and designed) as right chirality and/or double bond geometry or \"No isomers\" "
			sys.exit()
	
	if web_user_rightconfig_spec == "No" and web_user_delivery == "Pure":
		if chiral == True or d_bonds_out_ring > 0.0:
			print "ERROR| Your compound has unspecified chiral atom(s) and/or double bond(s), it can delivered as Mixture"
			sys.exit()

	if web_user_rightconfig_spec == "Yes":
		if chiral == False and d_bonds_out_ring == 0.0:
			print "ERROR| No chiral atoms or double bonds detected, please specify \"No isomers\" to the question on the right design or modify your compound"
			sys.exit()
		if chiral == True and specified_chirality == False:
			print "ERROR| Chiral atom(s) detected but the configuration has been not graphically specified"
			sys.exit()
		if d_bonds_out_ring > 0.0 and defined_EZ == False:
			print "ERROR| Double bond(s) detected but the geometry has been not graphically specified"
			sys.exit()

	return True

###################################################################################

import sys, time
sys.path.append( '/usr/programs/openbabel-2.4.1/lib64/python2.4/site-packages' )

import openbabel, pybel, pp

smart_def_file = '/home/www/chemotheca.unicz.it/var/pains_definitions.smarts'

ppservers = ()
job_server = pp.Server( 2, ppservers=ppservers )
# print "Starting pp with", job_server.get_ncpus(), "workers"


web = False

if len( sys.argv ) > 1:
	if sys.argv[1] == '-c':
		smile = sys.argv[2]
		complexity = IsComplex( smile )
		print complexity
		sys.exit()
	if sys.argv[1] == '-v':
		verbose_out = True
        	smile = sys.argv[2]
	elif sys.argv[1] == '-w':
		if len( sys.argv ) == 5:
			web = True
			verbose_out = False
			smile = sys.argv[2]
			web_user_rightconfig_spec = sys.argv[3]
			web_user_delivery = sys.argv[4]
		else:
			print "Web mod requires more command line arguments"
			sys.exit()
	else:
		verbose_out = False
		smile = sys.argv[1]
else:
	verbose_out = True
        smile = raw_input( 'Enter SMILE code: ' )


n_HBA              = 0
n_HBD              = 0
n_arom_atms        = 0
n_atms_in_rings    = 0
n_atms_connected   = 0
n_polar_H          = 0
n_non_polar_H      = 0
n_C                = 0
n_N                = 0
n_O                = 0
n_S                = 0
n_P                = 0
n_Cl               = 0
n_Br               = 0
n_I                = 0
n_F                = 0
n_carboxyl_O       = 0
n_phosphate_O      = 0
n_sulfate_O        = 0
n_nitro_O          = 0
n_amide_N          = 0
n_metals           = 0
n_aromatic_N_oxide = 0
n_primary_amide    = 0
n_secondary_amide  = 0
n_tertiary_amide   = 0
n_carboxyl_groups  = 0
n_phosphate_groups = 0
n_sulfate_groups   = 0
n_nitro_groups     = 0
n_ester_groups     = 0
n_carbonyl_groups  = 0
n_cheton_groups    = 0
n_aldehyde_groups  = 0
n_halogens         = 0
n_alcohol_O        = 0
n_pri_alcohol_O    = 0
n_sec_alcohol_O    = 0
n_ter_alcohol_O    = 0
n_ether_O          = 0
n_thioether_S      = 0
n_enol_O           = 0
n_phenol_O         = 0
n_not_chonsphalo   = 0
fpt                = ""
molsdf2D           = ""
molsdf3D           = ""
molmol2            = ""
pains_match        = "No"
pains_n_match      = 0
substr_name        = []
highlight          = []

obConversion = openbabel.OBConversion()
obConversion.SetInFormat("smi")

obmol = openbabel.OBMol()
obConversion.ReadString( obmol, smile )


pbmol = pybel.Molecule( obmol )

chiral, n_chiral_atms = ChiralAtoms( obmol) #
specified_chirality   = IsChiralitySpecified( smile, n_chiral_atms ) #
d_bonds       = pbmol.calcdesc( [ 'dbonds' ] ) [ 'dbonds' ]
d_bonds_in_ring, d_bonds_out_ring, defined_EZ = DoubleBondsAna( obmol, d_bonds, smile )

charge        = obmol.GetTotalCharge()
if charge != 0:
	salt = True
else:
	salt = False

if web:
	check_web_data( web_user_rightconfig_spec, web_user_delivery, chiral, specified_chirality, d_bonds_out_ring, defined_EZ )

obmol.DeleteHydrogens()
n_heavy_atoms = obmol.NumHvyAtoms()

formula       = obmol.GetFormula()
m_w           = obmol.GetMolWt()
n_bonds       = obmol.NumBonds()
s_bonds       = pbmol.calcdesc( [ 'sbonds' ] ) [ 'sbonds' ]
t_bonds       = pbmol.calcdesc( [ 'tbonds' ] ) [ 'tbonds' ]
a_bonds       = pbmol.calcdesc( [ 'abonds' ] ) [ 'abonds' ]
logP          = pbmol.calcdesc( [ 'logP' ] ) [ 'logP' ]   #
tpsa          = pbmol.calcdesc( [ 'TPSA' ] ) [ 'TPSA' ]   #
n_ha_rb       = obmol.NumRotors()

n_rings, n_aromatic_rings, n_aliphatic_rings = RingsCount( obmol )


aa_obmol      = obmol
aa_obmol.AddHydrogens()

n_all_atoms   = aa_obmol.NumAtoms()
n_aa_bonds    = aa_obmol.NumBonds()

n_bonds_to_H  = n_aa_bonds - n_bonds
n_hydrogens   = n_all_atoms - n_heavy_atoms


for obatom in openbabel.OBMolAtomIter( aa_obmol ):

	if charge == 0:
		n_connection = 0
		for obatm2 in openbabel.OBMolAtomIter( aa_obmol ):
			if obatom.IsConnected( obatm2 ):
				n_connection = n_connection + 1
		if n_connection > 0:
			n_atms_connected = n_atms_connected + 1

	if obatom.IsInRing():
		n_atms_in_rings = n_atms_in_rings + 1

	if obatom.IsPhosphorus():
		n_P = n_P + 1

	if obatom.IsSulfur():
		n_S = n_S + 1

	if obatom.IsNonPolarHydrogen():
		n_non_polar_H = n_non_polar_H + 1

	if obatom.IsPolarHydrogen():
		n_polar_H = n_polar_H + 1

	if obatom.IsAromatic():
		n_arom_atms = n_arom_atms + 1

	if obatom.IsHbondAcceptor():
		n_HBA = n_HBA + 1

	if obatom.IsHbondDonor():
		n_HBD = n_HBD + 1

	if obatom.IsCarbon():
		n_C = n_C + 1

	if obatom.IsNitrogen():
		n_N = n_N + 1

	if obatom.IsOxygen():
		n_O = n_O + 1

	if obatom.IsCarboxylOxygen():
		n_carboxyl_O = n_carboxyl_O + 1

	if obatom.IsPhosphateOxygen():
		n_phosphate_O = n_phosphate_O + 1

	if obatom.IsSulfateOxygen():
		n_sulfate_O = n_sulfate_O + 1

	if obatom.IsNitroOxygen():
		n_nitro_O = n_nitro_O + 1

	if obatom.IsAmideNitrogen():
		n_amide_N = n_amide_N + 1

	if obatom.IsAromaticNOxide():
		n_aromatic_N_oxide = n_aromatic_N_oxide + 1

	if obatom.IsMetal():
		n_metals = n_metals + 1
		if obatom.GetFormalCharge() == 0:
			print "ERROR| Metal atom(s) detected. Please specify its/their formal charge"
			sys.exit()

	if obatom.GetType() == 'Cl':
		n_Cl = n_Cl + 1 

	if obatom.GetType() == 'Br':
		n_Br = n_Br + 1 

	if obatom.GetType() == 'F':
		n_F = n_F + 1 

	if obatom.GetType() == 'I':
		n_I = n_I + 1 


if n_carboxyl_O > 0:
	n_carboxyl_groups = n_carboxyl_O / 2


if n_phosphate_O > 0:
	n_phosphate_groups = n_phosphate_O / 3


if n_sulfate_O > 0:
	n_sulfate_groups = n_sulfate_O / 3


if n_nitro_O > 0:
	n_nitro_groups = n_nitro_O / 2


n_halogens = n_Cl + n_Br + n_F + n_I


if n_O - ( n_carboxyl_O + n_phosphate_O + n_sulfate_O + n_nitro_O ) > 0:
	pri_alcohol_smarts = pybel.Smarts("[#6][CH2][OH]")
	sec_alcohol_smarts = pybel.Smarts("[#6][CH]([#6])[OH]")
	ter_alcohol_smarts = pybel.Smarts("[#6][C]([#6])([#6])[OH]")
	ether_smarts       = pybel.Smarts("[#6]O[#6]")
 	enol_smarts        = pybel.Smarts("[#6,#1]C([#6,#1])=C([#6,#1])[OH]")
	phenol_smarts      = pybel.Smarts("c[OH]")
	n_pri_alcohol_O    = len( pri_alcohol_smarts.findall( pbmol ) )
	n_sec_alcohol_O    = len( sec_alcohol_smarts.findall( pbmol ) )
	n_ter_alcohol_O    = len( ter_alcohol_smarts.findall( pbmol ) )
	n_ether_O          = len( ether_smarts.findall( pbmol ) )
	n_enol_O           = len( enol_smarts.findall( pbmol ) )
	n_phenol_O         = len( phenol_smarts.findall( pbmol ) )
	n_alcohol_O        = n_pri_alcohol_O + n_sec_alcohol_O + n_ter_alcohol_O


if n_S - n_sulfate_groups > 0:
	thioether_smarts = pybel.Smarts("[#6]S[#6]" )
	n_thioether_S = len( thioether_smarts.findall( pbmol ) )


for bond in openbabel.OBMolBondIter( obmol ):

	if bond.IsPrimaryAmide():
		n_primary_amide = n_primary_amide + 1

	if bond.IsSecondaryAmide():
		n_secondary_amide = n_secondary_amide + 1

	if bond.IsTertiaryAmide():
		n_tertiary_amide = n_tertiary_amide + 1

	if bond.IsEster() and my_IsNotCarboxylGroup( bond ):
		n_ester_groups = n_ester_groups + 1

#	if bond.IsCarbonyl() and my_IsNotCarboxylGroup( bond ):
#
#		if my_IsChetonGroup( bond ):
#			n_cheton_groups = n_cheton_groups + 1
#		
#		if my_IsAldehydeGroup( bond ):
#			n_aldehyde_groups = n_aldehyde_groups + 1
			

aldehyde_smarts = pybel.Smarts( "[CX3H1](=O)[#6]" )
n_aldehyde_groups = len( aldehyde_smarts.findall( pbmol ) )

cheton_smarts = pybel.Smarts( "[#6][CX3](=O)[#6]" )
n_cheton_groups = len( cheton_smarts.findall( pbmol ) )

n_carbonyl_groups = n_cheton_groups + n_aldehyde_groups


r3, r5         = Lipinski( n_HBA, n_HBD, logP, m_w ) #
logBB1         = logBB( logP, tpsa ) #


if logBB1 >= 0.3:
	CNS_Avail = True
else:
	CNS_Avail = False

if n_all_atoms != n_atms_connected:
	salt = True
else:
	salt = False

fpt = my_GetFingerPrint( obmol )

compute_molsdf2D = job_server.submit(my_2DBuilder, ("sdf", smile,), (), ())
molsdf2D = compute_molsdf2D()

compute_molsdf3D = job_server.submit(my_3DBuilder, ("sdf", smile,), (), ()) #sdf2smile, optimization
molsdf3D = compute_molsdf3D()

compute_molmol2  = job_server.submit(format_A_to_format_B, ("sdf", "mol2", molsdf3D,), (my_3DBuilder,), ())
molmol2 = compute_molmol2()

compute_inchi = job_server.submit(GetInChI, (smile,), (), ())
inchi = compute_inchi()

compute_inchikey = job_server.submit(GetInChIKey, (smile,), (), ())
inchikey = compute_inchikey()

# reconvert smiles due to strange string from JSME
compute_smile = job_server.submit(format_A_to_format_B, ("smi", "smi", smile,),(), ())
smile = compute_smile().rstrip()

# PAINS Check
pains_match, substr_name, highlight = PainsFinder( smile, smart_def_file )
pains_n_match = len(substr_name)

n_not_chonsphalomet = n_all_atoms - ( n_C + n_hydrogens + n_O + n_N + n_S + n_P + n_metals + n_halogens )


if verbose_out:
	print 'salt               :',salt
	print 'n_HBA              :',n_HBA
	print 'n_HBD              :',n_HBD
	print 'n_arom_atms        :',n_arom_atms
	print 'n_polar_H          :',n_polar_H
	print 'n_non_polar_H      :',n_non_polar_H
	print 'r3                 :',r3
	print 'r5                 :',r5
	print 'tpsa               :',tpsa
	print 'logBB1             :',logBB1
	print 'CNS_Avail          :',CNS_Avail
	print 'logP               :',logP
	print 'n_heavy_atoms      :',n_heavy_atoms  
	print 'n_all_atoms        :',n_all_atoms 
	print 'formula            :',formula 
	print 'm_w                :',m_w 
	print 'charge             :',charge 
	print 'chiral             :',chiral 
	print 'n_chiral_atms      :',n_chiral_atms 
	print 'specified_chirality:', specified_chirality
	print 'n_bonds    au      :',n_bonds,' aa:',n_aa_bonds 
	print 's_bonds    au      :',s_bonds
	print 'd_bonds    au      :',d_bonds
	print 'd_bonds_in_ring    :',d_bonds_in_ring
	print 'd_bonds_out_ring   :',d_bonds_out_ring
	print 'defined_EZ         :',defined_EZ
	print 't_bonds    au      :',t_bonds
	print 'a_bonds    au      :',a_bonds
	print 'n_heavy_atoms_rb   :',n_ha_rb
	print 'n_rings            :',n_rings 
	print 'InChI              :',inchi
	print 'InChIKey           :',inchikey
	print 'n_bonds_to_H       :',n_bonds_to_H
	print 'n_hydrogens        :',n_hydrogens
	print 'n_rings            :',n_rings
	print 'n_aliphatic_rings  :',n_aliphatic_rings
	print 'n_aromatic_rings   :',n_aromatic_rings
	print 'n_C                :',n_C 
	print 'n_N                :',n_N 
	print 'n_O                :',n_O 
	print 'n_S                :',n_S 
	print 'n_P                :',n_P 
	print 'n_Cl               :',n_Cl
	print 'n_Br               :',n_Br
	print 'n_F                :',n_F 
	print 'n_I                :',n_I 
	print 'n_enol_O           :',n_enol_O
	print 'n_phenol_O         :',n_phenol_O
	print 'n_alcohol_O        :',n_alcohol_O
	print 'n_pri_alcohol_O    :',n_pri_alcohol_O
	print 'n_sec_alcohol_O    :',n_sec_alcohol_O
	print 'n_ter_alcohol_O    :',n_ter_alcohol_O
	print 'n_ether_O          :',n_ether_O
	print 'n_thioether_S      :',n_thioether_S
	print 'n_carboxyl_O       :',n_carboxyl_O
	print 'n_carboxyl_groups  :',n_carboxyl_groups
	print 'n_phosphate_O      :',n_phosphate_O
	print 'n_phosphate_groups :',n_phosphate_groups
	print 'n_sulfate_O        :',n_sulfate_O 
	print 'n_sulfate_groups   :',n_sulfate_groups
	print 'n_nitro_O          :',n_nitro_O 
	print 'n_nitro_groups     :',n_nitro_groups
	print 'n_amide_N          :',n_amide_N 
	print 'n_aromatic_N_oxide :',n_aromatic_N_oxide 
	print 'n_metals           :',n_metals
	print 'n_primary_amide    :',n_primary_amide
	print 'n_secondary_amide  :',n_secondary_amide
	print 'n_tertiary_amide   :',n_tertiary_amide
	print 'n_ester_groups     :',n_ester_groups
	print 'n_carbonyl_groups  :',n_carbonyl_groups
	print 'n_cheton_groups    :',n_cheton_groups
	print 'n_aldehyde_groups  :',n_aldehyde_groups
	print 'n_halogens         :',n_halogens
	print 'fingerprint        :',fpt   
	print 'n_not_chonsphalomet:',n_not_chonsphalomet
	print 'SDF_3D             :',molsdf3D
	print 'SDF_2D             :',molsdf2D
	print 'MOL2               :',molmol2
	print 'SMILES             :',smile
	print 'PAINS              :',pains_match
	print 'PAINS-nMATCH       :',pains_n_match
	print 'PAINS-SUBSTRs      :',substr_name
	print 'ATOMS-MATCHING     :',highlight

else:

	substructures = ''
	atom_list     = ''
	if pains_n_match > 0:
		for i in range(0, pains_n_match):
			substructures = substructures + substr_name[i]
			for x in range(0, len(highlight[i])):
				atom_list = atom_list + str(highlight[i][x])
				if x < len(highlight[i]) - 1:
					atom_list = atom_list + ','
			if i < pains_n_match - 1:
				substructures = substructures + ';'
				atom_list = atom_list + ';'
	else:
		substructures = '0'
		atom_list     = '0'

  	print m_w,"|",tpsa,"|",formula,"|",r3,"|",r5,"|",logP,"|",n_HBA,"|",n_HBD,"|",logBB1,"|",CNS_Avail,"|",n_heavy_atoms,"|",n_all_atoms,"|",n_bonds,"|",n_aa_bonds,"|",n_arom_atms,"|",n_polar_H,"|",n_non_polar_H,"|",n_all_atoms,"|",charge,"|",chiral,"|",salt,"|",n_chiral_atms,"|",specified_chirality,"|",s_bonds,"|",d_bonds,"|",d_bonds_in_ring,"|",d_bonds_out_ring,"|",defined_EZ,"|",t_bonds,"|",a_bonds,"|",n_ha_rb,"|",n_rings,"|",inchi,"|",inchikey,"|",n_bonds_to_H,"|",n_hydrogens,"|",n_rings,"|",n_aliphatic_rings,"|",n_aromatic_rings,"|",n_C,"|",n_N,"|",n_O,"|",n_S,"|",n_P,"|",n_Cl,"|",n_Br,"|",n_F,"|",n_I,"|",n_carboxyl_O,"|",n_carboxyl_groups,"|",n_phosphate_O,"|",n_phosphate_groups,"|",n_sulfate_O,"|",n_sulfate_groups,"|",n_nitro_O,"|",n_nitro_groups,"|",n_amide_N,"|",n_aromatic_N_oxide,"|",n_metals,"|",n_primary_amide,"|",n_secondary_amide,"|",n_tertiary_amide,"|",n_ester_groups,"|",n_carbonyl_groups,"|",n_cheton_groups,"|",n_aldehyde_groups,"|",n_halogens,"|",fpt,"|",molsdf3D,"|",n_alcohol_O,"|",n_pri_alcohol_O,"|",n_sec_alcohol_O,"|",n_ter_alcohol_O,"|",n_ether_O,"|",n_thioether_S,"|",n_enol_O,"|",n_phenol_O,"|",n_not_chonsphalomet,"|",molsdf2D,"|",molmol2,"|",smile,"|",pains_match,"|",pains_n_match,"|",substructures,"|",atom_list
