#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compute the adsorbed geometry of a molecule on a surface using meaningful geometrical constraints"""

# Importing librairies

from __future__ import print_function
import sys
import argparse
import numpy as np
import ase
from ase import io
from ase.io import read
from ase.io import write


# Module metadata

__author__ = "Ruben Staub"
__copyright__ = "Copyright © Laboratoire de Chimie, ENS de Lyon"
__credits__ = ["Ruben Staub", "Sarah Blanck", "Danish Kaur Pannu", "Carine Seraphim"]
__license__ = "GNU LGPLv3"
__version__ = "1.0"
__maintainer__ = "Ruben Staub"
__email__ = "ruben.staub@ens-lyon.fr"
__status__ = "Release"


# Functions definitions

def get_proper_angle(v1, v2, ref=None, degrees=True):
	"""Computes the signed angle between v2 and v1, counterclockwise when facing ref.
	
	Parameters:
		v1 (numpy.ndarray): First vector (used as reference).
		
		v2 (numpy.ndarray): Second vector.
		
		ref (numpy.ndarray) (optional): Orthogonal vector to both v1 and v2, along which the sign of the rotation is defined (i.e. positive if counterclockwise angle when facing ref).
			Note: If None, angle will be given between 0 and 180°.
		
		degrees (boolean) (optional): Control whether the angle is returned in degrees (=True by default) or in radians.
	
	Example:
		In a direct basis, get_proper_angle([1,0,0], [0,1,0], [0,0,1])
	
	Note:
		This function is required for robust computations of angles, since np.arccos handle very poorly round-off errors natively...
		Example: ase.geometry.get_angles([[0,3*np.pi,3*np.pi]],[[0,np.pi,np.pi]]) might just fail...
	"""
#	print(v1, v2, ase.geometry.get_angles([v1],[v2]))
	# Computes dot product of normalized vectors
	norm_dot = np.dot(v1/np.linalg.norm(v1), v2/np.linalg.norm(v2))
	
	# Computes angle
	angle = np.arccos(np.sign(norm_dot) if np.abs(norm_dot) >= 1 else norm_dot)
	
	if ref is not None:
		# Give sign according to ref direction
		angle *= (1 if np.dot(np.cross(v1, v2), ref) >= 0 else -1) # Trick to get the sign (since the angle is initially found between 0 and 180°)
	
	return(angle*180/np.pi if degrees else angle)

def transform_adsorbate(molecule, surface, atom1_mol, atom2_mol, atom3_mol, atom1_surf, atom2_surf, bond_vector, bond_angle_target, dihedral_angle_target=None, mol_dihedral_angle_target=None):
	"""Performs translation and rotation of an adsorbate defined by an adsorption bond length, direction, angle and dihedral angle
	
	Parameters:
		molecule (ase.Atoms): The molecule to adsorbate.
		
		surface (ase.Atoms): The surface onto which the molecule must be adsorbed.
		
		atom1_mol (ase.Atom): The atom of the adsorbate that will be bound to the surface.
		
		atom2_mol (ase.Atom): An other atom of the adsorbate used to define the adsorption bond angle, and the dihedral adsorption angle.
		
		atom3_mol (ase.Atom): An third atom of the adsorbate used to define the adsorbate dihedral angle.
		
		atom1_surf (ase.Atom): The atom of the surface that will be bound to the molecule.
		
		atom2_surf (ase.Atom): An other atom of the surface used to define the dihedral adsorption angle.
		
		bond_vector (numpy.ndarray): The adsorption bond desired.
			Details: offset = vect(atom1_surf, atom1_mol)
		
		bond_angle_target (float or int): The adsorption bond angle desired (in degrees).
			Details: bond_angle_target = angle(atom1_surf, atom1_mol, atom2_mol)
		
		dihedral_angle_target (float or int): The dihedral adsorption angle desired (in degrees).
			Details: dihedral_angle_target = dihedral(atom2_surf, atom1_surf, atom1_mol, atom2_mol)
				The rotation vector is facing the adsorbate from the surface (i.e. counterclockwise rotation when facing the surface (i.e. view from top))
		
		mol_dihedral_angle_target (float or int): The adsorbate dihedral angle desired (in degrees).
			Details: mol_dihedral_angle_target = dihedral(atom1_surf, atom1_mol, atom2_mol, atom3_mol)
				The rotation vector is facing atom2_mol from atom1_mol
	
	Returns:
		None (the `molecule` object is modified in-place)
	"""
	# Retrieve bonds of interest
	bond_surf = atom2_surf.position - atom1_surf.position
	bond_inter = atom1_mol.position - atom1_surf.position
	bond_mol = atom2_mol.position - atom1_mol.position
	bond2_mol = atom3_mol.position - atom2_mol.position
	
	# Check if dihedral angles can be defined
	do_dihedral = dihedral_angle_target is not None
	do_mol_dihedral = mol_dihedral_angle_target is not None
	dihedral_use_mol2 = False
	# Check if bond_surf and bond_inter are not aligned
	if np.allclose(np.cross(bond_surf, bond_inter), 0):
		print("Warning: Surface atoms are incompatible with adsorption direction/bond. An adsorption dihedral angle cannot be defined.", file=sys.stderr)
		do_dihedral = False
	# Check if requested bond angle is not flat
	if np.isclose((bond_angle_target + 90)%180 - 90, 0):
		print("Warning: Requested bond angle is flat. Only a single dihedral angle can be defined (atom2_surf, atom1_surf, atom2_mol, atom3_mol).", file=sys.stderr)
		do_mol_dihedral = False
		dihedral_use_mol2 = True
		print("Warning: Dihedral adsorption angle rotation will be perfomed with ({}, {}, {}, {}).".format(atom2_surf.index, atom1_surf.index, atom2_mol.index, atom3_mol.index), file=sys.stderr)
	# Check if bond_mol and bond2_mol are not aligned
	if np.allclose(np.cross(bond_mol, bond2_mol), 0):
		print("Warning: Adsorbates atoms are aligned. An adsorbate dihedral angle cannot be defined.", file=sys.stderr)
		do_mol_dihedral = False
	if not do_dihedral:
		print("Warning: Dihedral adsorption angle rotation will not be performed.", file=sys.stderr)
	if not do_mol_dihedral:
		print("Warning: Adsorbate dihedral angle rotation will not be performed.", file=sys.stderr)
	
	###########################
	#       Translation       #
	###########################
	
	# Compute and apply translation of adsorbate
	translation = bond_vector - bond_inter
	molecule.translate(translation)
	
	# Update adsorption bond
	bond_inter = atom1_mol.position - atom1_surf.position
	
	# Check if translation was successful
	if np.allclose(bond_inter, bond_vector):
		print("Translation successfully applied (error: ~ {:.5g} unit length)".format(np.linalg.norm(bond_inter-bond_vector)))
	else:
		raise AssertionError('An unknown error occured during the translation')
	
	###########################
	#   Bond angle rotation   #
	###########################
	
	# Compute rotation vector
	rotation_vector = np.cross(-bond_inter, bond_mol)
	if np.allclose(rotation_vector, 0, atol=1e-5):
		# If molecular bonds are aligned, any vector orthogonal to bond_inter can be used
		# Such vector can be found as the orthogonal rejection of either X-axis, Y-axis or Z-axis onto bond_inter (since they cannot be all aligned)
		non_aligned_vector = np.zeros(3)
		non_aligned_vector[np.argmin(np.abs(bond_inter))] = 1 # Select the most orthogonal axis (lowest dot product)
		rotation_vector = non_aligned_vector - np.dot(non_aligned_vector, bond_inter)/np.dot(bond_inter, bond_inter) * bond_inter
	
	# Retrieve current bond angle
	bond_angle_ini = get_proper_angle(-bond_inter, bond_mol, rotation_vector)
	
	# Apply rotation to reach desired bond_angle
	molecule.rotate(bond_angle_target-bond_angle_ini, v=rotation_vector, center=atom1_mol.position)
	
	# Update molecular bonds
	bond_mol = atom2_mol.position - atom1_mol.position
	bond2_mol = atom3_mol.position - atom2_mol.position
	
	# Check if rotation was successful
	bond_angle = get_proper_angle(-bond_inter, bond_mol)
	if np.isclose((bond_angle - bond_angle_target + 90)%180 - 90, 0, atol=1e-3) and np.allclose(atom1_mol.position - atom1_surf.position, bond_inter):
		print("Rotation successfully applied (error: {:.5f}°)".format((bond_angle - bond_angle_target + 90)%180 - 90))
	else:
		raise AssertionError('An unknown error occured during the rotation')
	
	###########################
	# Dihedral angle rotation #
	###########################
	
	# Perform dihedral rotation if possible
	if do_dihedral:
		# Retrieve current dihedral angle (by computing the angle between the orthogonal rejection of bond_surf and bond_mol onto bond_inter)
		bond_inter_inner = np.dot(bond_inter, bond_inter)
		bond_surf_reject = bond_surf - np.dot(bond_surf, bond_inter)/bond_inter_inner * bond_inter
		if dihedral_use_mol2:
			bond_mol_reject = bond2_mol - np.dot(bond2_mol, bond_inter)/bond_inter_inner * bond_inter
		else:
			bond_mol_reject = bond_mol - np.dot(bond_mol, bond_inter)/bond_inter_inner * bond_inter
		dihedral_angle_ini = get_proper_angle(bond_surf_reject, bond_mol_reject, bond_inter)
		
		# Apply dihedral rotation along bond_inter
		molecule.rotate(dihedral_angle_target-dihedral_angle_ini, v=bond_inter, center=atom1_mol.position)
		
		# Update molecular bonds
		bond_mol = atom2_mol.position - atom1_mol.position
		bond2_mol = atom3_mol.position - atom2_mol.position
		
		# Check if rotation was successful
		# Check dihedral rotation
		if dihedral_use_mol2:
			bond_mol_reject = bond2_mol - np.dot(bond2_mol, bond_inter)/bond_inter_inner * bond_inter
		else:
			bond_mol_reject = bond_mol - np.dot(bond_mol, bond_inter)/bond_inter_inner * bond_inter
		dihedral_angle = get_proper_angle(bond_surf_reject, bond_mol_reject, bond_inter)
		# Check bond rotation is unmodified
		bond_angle = get_proper_angle(-bond_inter, bond_mol)
		if np.isclose((dihedral_angle - dihedral_angle_target + 90)%180 - 90, 0, atol=1e-3) and np.isclose((bond_angle - bond_angle_target + 90)%180 - 90, 0, atol=1e-5) and np.allclose(atom1_mol.position - atom1_surf.position, bond_inter):
			print("Dihedral rotation successfully applied (error: {:.5f}°)".format((dihedral_angle - dihedral_angle_target + 90)%180 - 90))
		else:
			raise AssertionError('An unknown error occured during the dihedral rotation')
	
	#####################################
	# Adsorbate dihedral angle rotation #
	#####################################
	
	# Perform adsorbate dihedral rotation if possible
	if do_mol_dihedral:
		# Retrieve current adsorbate dihedral angle (by computing the angle between the orthogonal rejection of bond_inter and bond2_mol onto bond_mol)
		bond_mol_inner = np.dot(bond_mol, bond_mol)
		bond_inter_reject = -bond_inter - np.dot(-bond_inter, bond_mol)/bond_mol_inner * bond_mol
		bond2_mol_reject = bond2_mol - np.dot(bond2_mol, bond_mol)/bond_mol_inner * bond_mol
		dihedral_angle_ini = get_proper_angle(bond_inter_reject, bond2_mol_reject, bond_mol)
		
		# Apply dihedral rotation along bond_mol
		molecule.rotate(mol_dihedral_angle_target-dihedral_angle_ini, v=bond_mol, center=atom1_mol.position)
		
		# Update molecular bonds
		bond_mol = atom2_mol.position - atom1_mol.position
		bond2_mol = atom3_mol.position - atom2_mol.position
		
		# Check if rotation was successful
		# Check adsorbate dihedral rotation
		bond_mol_inner = np.dot(bond_mol, bond_mol)
		bond2_mol_reject = bond2_mol - np.dot(bond2_mol, bond_mol)/bond_mol_inner * bond_mol
		mol_dihedral_angle = get_proper_angle(bond_inter_reject, bond2_mol_reject, bond_mol)
		# Check dihedral rotation
		bond_mol_reject = bond_mol - np.dot(bond_mol, bond_inter)/bond_inter_inner * bond_inter
		dihedral_angle = get_proper_angle(bond_surf_reject, bond_mol_reject, bond_inter)
		# Check bond rotation is unmodified
		bond_angle = get_proper_angle(-bond_inter, bond_mol)
		if np.isclose((mol_dihedral_angle - mol_dihedral_angle_target + 90)%180 - 90, 0, atol=1e-3) and np.isclose((dihedral_angle - dihedral_angle_target + 90)%180 - 90, 0, atol=1e-5) and np.isclose((bond_angle - bond_angle_target + 90)%180 - 90, 0, atol=1e-5) and np.allclose(atom1_mol.position - atom1_surf.position, bond_inter):
			print("Adsorbate dihedral rotation successfully applied (error: {:.5f}°)".format((mol_dihedral_angle - mol_dihedral_angle_target + 90)%180 - 90))
		else:
			raise AssertionError('An unknown error occured during the adsorbate dihedral rotation')


# Main program

if __name__ == '__main__':
	# Read from command line arguments
	parser = argparse.ArgumentParser()
	parser.add_argument("base_input_filename", help="filename of base input")
	parser.add_argument("new_input_filename", help="filename of new input")
	
	parser.add_argument("base1_index", help="index of base1 atom", type=int)
	parser.add_argument("base2_index", help="index of base2 atom", type=int)
	parser.add_argument("new1_index", help="index of new1 atom", type=int)
	parser.add_argument("new2_index", help="index of new2 atom", type=int)
	parser.add_argument("new3_index", help="index of new3 atom", type=int)
	
	parser.add_argument("-bx", "--bond_vector_x", help="X-component of requested bond vector", type=float, default=0.0)
	parser.add_argument("-by", "--bond_vector_y", help="Y-component of requested bond vector", type=float, default=0.0)
	parser.add_argument("bond_vector_z", help="Z-component of requested bond vector", type=float)
	
	parser.add_argument("bond_angle", help="requested bond angle (base1, new1, new2)", type=float)
	parser.add_argument("bond_dihedral_angle", help="requested bond dihedral angle (base2, base1, new1, new2)", type=float)
	parser.add_argument("new_dihedral_angle", help="requested new dihedral angle (base1, new1, new2, new3)", type=float)
	
	parser.add_argument("output_filename", help="filename of concatenated output")
	
	parser.add_argument("--force_cartesian", help="force writing cartesian coordinates (even if PBC is detected)", action="store_true")
	parser.add_argument("--vasp5", help="use VASP5 file format", action="store_true")
	parser.add_argument("--sort", help="sort atomic coordinates (useful for VASP format)", action="store_true")
	
	args = parser.parse_args()
	
	
	# Retrieve objects
	base_structure = read(args.base_input_filename)
	print("Base structure ({} atoms) successfully retrieved from {}".format(base_structure.get_number_of_atoms(), args.base_input_filename))
	new_structure = read(args.new_input_filename)
	print("New structure ({} atoms) successfully retrieved from {}\n".format(new_structure.get_number_of_atoms(), args.new_input_filename))
	
	base1_atom = base_structure[args.base1_index]
	print("Selected base1 atom is {}".format(base1_atom))
	base2_atom = base_structure[args.base2_index]
	print("Selected base2 atom is {}".format(base2_atom))
	new1_atom = new_structure[args.new1_index]
	print("Selected new1 atom is {}".format(new1_atom))
	new2_atom = new_structure[args.new2_index]
	print("Selected new2 atom is {}".format(new2_atom))
	new3_atom = new_structure[args.new3_index]
	print("Selected new3 atom is {}\n".format(new3_atom))
	
	bond_vector = np.array([args.bond_vector_x, args.bond_vector_y, args.bond_vector_z])
	print("Requested bond vector (base1 -> new1) is {}".format(bond_vector))
	print("Requested bond angle (base1, new1, new2) is {}".format(args.bond_angle))
	print("Requested bond dihedral angle (base2, base1, new1, new2) is {}".format(args.bond_dihedral_angle))
	print("Requested new dihedral angle (base1, new1, new2, new3) is {}\n".format(args.new_dihedral_angle))
	
	write_options = dict()
	if args.vasp5:
		write_options['vasp5'] = True
	if args.sort:
		write_options['sort'] = True
	if args.force_cartesian:
		write_options['direct'] = False
	print("Additional options requested for ase.io.write: {}".format(write_options))
	
	
	# Perform transformation of new structure geometry
	print("\nPerfoming transformations of the new structure to reach desired properties:")
	transform_adsorbate(new_structure, base_structure, new1_atom, new2_atom, new3_atom, base1_atom, base2_atom, bond_vector, args.bond_angle, dihedral_angle_target=args.bond_dihedral_angle, mol_dihedral_angle_target=args.new_dihedral_angle)
	
	# Concatenate structures
	final = base_structure + new_structure
	
	# Check if periodic boundary conditions (PBC) is on
	pbc_bool = final.get_pbc().any()
	# Use direct coordinates output if PBC are on (unless prevented by user)
	# This is due to a bug of Molden: POSCAR files are written in cartesian by default, but the "Cartesian" keyword does not seem well handled by Molden...
	if pbc_bool:
		write_options['direct'] = write_options.get('direct', True)
		if write_options['direct']:
			print("\nPeriodic boundaries conditions detected, output will be written in direct coordinates (due to Molden compatibilities)")
		else:
			print("\nPeriodic boundaries conditions detected, cartesian coordinates requested")
			print("Warning: PBC is detected but cartesian coordinates are requested, expecting Molden issues for POSCAR output files...", file=sys.stderr)
		if not final.get_pbc().all():
			print("Warning: anisotropic PBC is not yet supported by ASE, expecting uncoherent geometries...", file=sys.stderr)
	
	# Write combined geometry
	write(args.output_filename, final, **write_options)
	print("\nCombined geometry (base+new) successfully written into '{}'".format(args.output_filename))
