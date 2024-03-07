#! /usr/bin/env python
# 2021 -- 2023 - Siebe Vanlommel
from molmod.io import dump_chk, load_chk
import molmod.graphs
from molmod.molecules import Molecule
from molmod.molecular_graphs import MolecularGraph
from molmod.graphs import GraphSearch, CustomPattern
from molmod.molecular_graphs import HasAtomNumber, HasNumNeighbors, HasNeighborNumbers, HasNeighbors
from molmod.graphs import CritAnd, CritNot, CritOr
from molmod.periodic import periodic as pt
from molmod.units import angstrom
from yaff import System as YaffSystem
from yaff.log import log
import sys, os, copy, numpy as np
from collections import Counter
import itertools
from ase.io import read

log.set_level(log.silent)

fn_xyz = read('system.extxyz',format='extxyz')
cell = fn_xyz.get_cell()

#print(cell)
sys = YaffSystem.from_file('system.xyz',rvecs = cell*angstrom)
global_nums = sys.numbers
sys.to_file('system.chk')

from molmod.unit_cells import UnitCell

def get_pattern(target, cof, a_count = 0):
	targ = Molecule.from_file(target)
	targ_nums = targ.numbers
	targ_graph = MolecularGraph.from_geometry(targ)
	cof = YaffSystem.from_file(cof, rvecs = np.array(cell)*angstrom)
	cof.detect_bonds()
	cof_nums = cof.numbers
	cof_graph = MolecularGraph(cof.bonds, cof.numbers)
	pattern = CustomPattern(targ_graph)
	gs = GraphSearch(pattern)

	count = 0
	out = {}
	for match in gs(cof_graph):
		count += 1
		out[count] = match.forward

	map = {}
	unique = []
	import quickff
	import quickff.tools as tls

	tls.set_ffatypes(cof,'high')
	#print(out)
	labels = {} # dict where key=(index within pattern) and val=(list of cof indices belonging to that class)
	for ii in range(len(targ.numbers)):
		labels[ii] = []
	for c in out.keys():
		map[c] = [out[c][a] for a in out[c].keys()]
		if not set(map[c]) in unique:
			unique.append(set(map[c]))
			# now assign cof indices to pattern indices
			for a in out[c].keys():
				labels[a].append(out[c][a])

	#print(unique)

	for m in unique:
		a_count += len(m)

	return a_count, unique, labels, targ_nums, cof


if __name__ == '__main__':
	# this yields 6 nodes as defined in TAPD_reduced.xyz
	atomcount,nodes,node_labels,node_nums,cof_sys = get_pattern('node.xyz', 'system.xyz')
	print('number of nodes: ', len(nodes))
	# this yields 12 linkers as defined in node.xyz
	#print('Warning: 2-coordinate nitrogen are double accounted for (in both TAPD and Me)')
	atomcount,links,linker_labels,linker_nums,cof_sys = get_pattern('linker.xyz', 'system.xyz', atomcount)
	print('number of linkers: ',len(links))
	#print('atom count = ',atomcount,', which is due to double counting of N')
	print(' ')
	print('classification')
	print('  ')
	print('node index        corresponding cof indices')
	grouped = {}
	for n_i,cof_is in node_labels.items():
		if node_nums[n_i] == 6:
			assert np.all([global_nums[cof_i] == 6 for cof_i in cof_is])
			#print(n_i,'       ', *[(findex,cof_sys.ffatypes[cof_sys.ffatype_ids[findex]]) for findex in cof_is])
			grouped['node-{}'.format(n_i)] = [[findex,cof_sys.ffatypes[cof_sys.ffatype_ids[findex]]] for findex in cof_is]

	print(' ')
	print('linker index      corresponding (cof indx,ffatype)')
	for l_i,cof_is in linker_labels.items():
		if linker_nums[l_i] == 6:
			assert np.all([global_nums[cof_i] == 6 for cof_i in cof_is])
			#print(l_i,'       ', *[(findex,cof_sys.ffatypes[cof_sys.ffatype_ids[findex]]) for findex in cof_is])
			grouped['linker-{}'.format(l_i)] = [[findex,cof_sys.ffatypes[cof_sys.ffatype_ids[findex]]] for findex in cof_is]
	writefile = open('linkers_classified.dat', 'w')
	writefile.close()
	writefile = open('nodes_classified.dat', 'w')
	writefile.close()
	#print('ffatypes in cof: ',cof_sys.ffatypes)
	print(grouped)
	for label, sites in grouped.items():
		kind, number= label.split('-')
		if kind == 'node': filename = 'nodes'
		else: filename = 'linkers'
		writefile = open('{}_classified.dat'.format(filename), 'a')
		group_save = np.array([number,[site[0] for site in sites]]).flatten()
		towrite = []
		towrite.append(str(group_save[0]))
		for item in group_save[1]:
			towrite.append(str(item))
		writefile.write(' '.join(towrite) + '\n')
		writefile.close()
	# save
	#import pickle
	#dumpfile = open("carbons_classified.dic", 'wb')
	#pickle.dump(grouped,dumpfile)
	#dumpfile.close()
