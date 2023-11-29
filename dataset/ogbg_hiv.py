from ogb.graphproppred import GraphPropPredDataset
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit import Chem
import numpy as np
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import Image
import json

def edge_index_to_mol(node_features, edge_index):
    mol = Chem.RWMol()
    for atomic_num in node_features:
        mol.AddAtom(Chem.Atom(int(atomic_num)))
    for src, dst in zip(edge_index[0], edge_index[1]):
        if not mol.GetBondBetweenAtoms(int(src), int(dst)):
            mol.AddBond(int(src), int(dst), Chem.BondType.SINGLE)
    mol = mol.GetMol()
    return mol

def count_six_membered_rings(mol):
    Chem.SanitizeMol(mol)
    num_membered_rings=defaultdict(int)
    ring_info = mol.GetRingInfo()
    for ring in ring_info.AtomRings():
        num_membered_rings[len(ring)]+=1
    return num_membered_rings

def draw_mol_to_svg(mol):
    d = Draw.MolDraw2DSVG(224, 224)
    opts = d.drawOptions()
    for i in range(mol.GetNumAtoms()):
        opts.atomLabels[i] = ""
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = "<svg xmlns='http://www.w3.org/2000/svg' width='224px' height='224px' viewBox='0 0 224 224'>\n"
    svg += "<rect style='fill:white' width='224.0' height='224.0' x='0.0' y='0.0'> </rect>\n"
    for path in d.GetDrawingText().split('class')[1:]:
        svg+= "<path d='M"+path.split("' d='M")[1].split("' style")[0]+"' style='stroke:black' />\n"
    svg+="</svg>"
    return svg
  
dataset = GraphPropPredDataset(name='ogbg-molhiv')# Select a graph from the dataset
for idx,graph in enumerate(dataset):
    edge_index = graph[0]['edge_index']
    node_features = np.zeros(graph[0]['num_nodes'])
    mol = edge_index_to_mol(node_features, edge_index)
    rings = count_six_membered_rings(mol)
    if len(rings.keys())==1 and 6 in rings.keys():
        output = {
            "idx" : idx,
            "edge_index" : ', '.join(f'({i}, {j})' for (i,j) in edge_index.T.tolist()),
            "svg" : draw_mol_to_svg(mol),
            "6 cycle" : rings[6]
        }
        with open(f'eval/ogbg-hiv_{rings[6]}.jsonl', "a") as jsonl_file:
            jsonl_file.write(json.dumps(output) + "\n")
