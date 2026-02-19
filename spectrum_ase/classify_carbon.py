# Siebe Vanlommel 2026--..
# molmod+yaff-free version partly made by Gemini 3
import networkx as nx
from ase.io import read
from ase.neighborlist import natural_cutoffs, NeighborList

def atoms_to_graph(atoms, cutoff_mult=1.1):
    """Converts ASE atoms to a NetworkX graph using covalent radii."""
    cutoffs = natural_cutoffs(atoms, mult=cutoff_mult)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        G.add_node(i, element=atom.symbol)
        
    for i in range(len(atoms)):
        indices, _ = nl.get_neighbors(i)
        for n in indices:
            G.add_edge(i, n)
    return G

def get_unique_assignments(G_sys, G_block, assigned_set):
    """
    Finds subgraph matches and assigns carbon atoms.
    Atoms in the block that are already in assigned_set are skipped,
    allowing linkers to share interface atoms with nodes.
    """
    node_match = lambda n1, n2: n1['element'] == n2['element']
    GM = nx.algorithms.isomorphism.GraphMatcher(G_sys, G_block, node_match=node_match)
    
    # Store results: {block_atom_index: [list of system_indices]}
    mappings = {i: [] for i, data in G_block.nodes(data=True) if data['element'] == 'C'}
    
    for match in GM.subgraph_isomorphisms_iter():
        # match: {sys_idx: block_idx}
        for s_idx, b_idx in match.items():
            if b_idx in mappings and s_idx not in assigned_set:
                mappings[b_idx].append(s_idx)
                assigned_set.add(s_idx)  # claim the atom

    return mappings


if __name__ == "__main__":
    # Load files
    system = read("system.xyz")
    node_ref = read("node.xyz")
    linker_ref = read("linker.xyz")

    # Build graphs
    print("Building graphs (this may take a moment for large systems)...")
    G_sys = atoms_to_graph(system)
    G_node = atoms_to_graph(node_ref)
    G_link = atoms_to_graph(linker_ref)

    # This set tracks atoms that have been 'claimed'
    global_assigned_atoms = set()

    # 1. Process Nodes first (usually the more rigid part of the framework)
    print("Fitting Nodes...")
    node_results = get_unique_assignments(G_sys, G_node, global_assigned_atoms)

    # 2. Process Linkers from the remaining atoms
    print("Fitting Linkers...")
    link_results = get_unique_assignments(G_sys, G_link, global_assigned_atoms)

    # Save results
    def save_classified(filename, data):
        with open(filename, 'w') as f:
            for b_idx in sorted(data.keys()):
                indices = [str(int(i)) for i in sorted(data[b_idx])]
                f.write(f"{b_idx} {' '.join(indices)}\n")

    save_classified("nodes_classified.dat", node_results)
    save_classified("linkers_classified.dat", link_results)

    # Verification
    total_sys_c = len([a for a in system if a.symbol == 'C'])
    n_c = sum(len(v) for v in node_results.values())
    l_c = sum(len(v) for v in link_results.values())
    
    print(f"\nSuccess!")
    print(f"Node Carbon Count: {n_c}")
    print(f"Linker Carbon Count: {l_c}")
    print(f"Total classified: {n_c + l_c} / {total_sys_c}")
    
    # Check for equality in groups
    if node_results:
        counts = [len(v) for v in node_results.values()]
        print(f"Node indices per group: {min(counts)} to {max(counts)} (Should be equal)")
    if link_results:
        counts = [len(v) for v in link_results.values()]
        print(f"Linker indices per group: {min(counts)} to {max(counts)} (Should be equal)")
