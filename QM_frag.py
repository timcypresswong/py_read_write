import argparse
import MDAnalysis as mda
from rdkit import Chem
from openbabel import openbabel
import networkx as nx
import pandas as pd
from rdkit.Chem import AllChem
from rdkit.Chem import rdEHTTools
import networkx as nx
import math
from openmm.app import PDBFile, ForceField, Simulation, Modeller
from openmm import CustomExternalForce, Platform, VerletIntegrator, LangevinIntegrator
from openmm import unit as openmm_unit
import openmm
from openmm import Vec3
from openff.toolkit import Molecule as OFFMolecule

def parser_args():
	parser = argparse.ArgumentParser(description='This is the select the molecular fragment based on selection method, this will complement the selected fragments to the C-C bond')
	parser.add_argument('filename', help='input pdb name')
	parser.add_argument('selection', help='specific selection command with quatation mark, i.e. \"(around 5 resname LIG) or (resname LIG)\" or \"same segid as (around 5 resname LIG)\", etc.')
	parser.add_argument('--sep', help='define the seperation mark for output, default is space', type=str)
	parser.add_argument('--base', help='determine whether the selected atoms indices are 0 or 1-based, default is 0-based', type=int)
	parser.add_argument('--output_pdb', help='specify the outputname of the pdb', type=str)
	args = parser.parse_args()
	return args


def judge_charge_multiplicity(num_electron, partial_charge, multiplicity = 1):
    floorcharge = math.floor(partial_charge)
    ceilcharge = math.ceil(partial_charge)
    num_unpair_electron = multiplicity - 1
    if ( - floorcharge + num_electron - num_unpair_electron ) % 2 == 0:
        return floorcharge
    else:
        return ceilcharge


def calc_Gasteiger_Charge(RDmol):
    charge_list = []
    AllChem.ComputeGasteigerCharges(RDmol, nIter = 24)
    for atom in RDmol.GetAtoms():
        g_charge = atom.GetProp('_GasteigerCharge')
        # 氢原子总电荷（含隐式氢贡献）[reference:6]
        # h_charge = atom.GetProp('_GasteigerHCharge')
        # print(f'原子 {atom.GetIdx()} 的 Gasteiger 电荷: {g_charge}')
        charge_list.append(g_charge)

    return charge_list

def calc_Formal_Charge(RDmol):
    charge_list = []
    for atom in RDmol.GetAtoms():
        g_charge = atom.GetFormalCharge()
        charge_list.append(g_charge)
    return charge_list


def calc_EHT_Charge(RDmol):
    charge_list = []
    _, res = rdEHTTools.RunMol(RDmol)
    charge_list = res.GetAtomicCharges()
    return charge_list


def get_Atomic_Num(RDmol):
    atomic_num_list = []
    for atom in RDmol.GetAtoms():
        atomic_num_list.append( atom.GetAtomicNum() )
    return atomic_num_list


def Merge_two_graphs_keep_high_BO(G1, G2):
    G = nx.Graph()
    edges = []
    all_edges, w1, w2 = aligned_edge_weights(G1, G2)
    for e, d1, d2 in zip( all_edges, w1, w2 ):

        d = max(d1, d2)
        edges.append(  tuple( [e[0], e[1], d]  )  )
    G.add_weighted_edges_from(edges)
    return G

def RDMol_to_Graph(RD_mol):
    G = nx.Graph()
    # 添加带权重的边
    edges = []
    for bond in RD_mol.GetBonds():
        bond_idx = bond.GetIdx()
        begin_atom_idx = bond.GetBeginAtomIdx()
        end_atom_idx = bond.GetEndAtomIdx()
        bond_order = int(bond.GetBondTypeAsDouble())
        edges.append( tuple([begin_atom_idx, end_atom_idx, bond_order]) )
    G.add_weighted_edges_from(edges)

    for atom in RD_mol.GetAtoms():
        if atom.GetTotalValence() == 0 or atom.GetDegree() == 0:
            G.add_node(atom.GetIdx())
    return G

def OBMol_to_Graph(OB_mol):
    G = nx.Graph()
    # 添加带权重的边
    edges = []
    for bond in openbabel.OBMolBondIter(OB_mol):
        bond_idx = bond.GetIdx()
        begin_atom_idx = bond.GetBeginAtom().GetIdx()
        end_atom_idx = bond.GetEndAtom().GetIdx()
        bond_order = bond.GetBondOrder()
        edges.append( tuple([begin_atom_idx - 1, end_atom_idx - 1, bond_order]) )

    G.add_weighted_edges_from(edges)


    for atom in openbabel.OBMolAtomIter(OB_mol):
        if atom.GetExplicitValence() == 0 or atom.GetExplicitDegree() == 0:
            G.add_node(atom.GetIdx() - 1)

    return G

def get_weight(G, edge, weight_key = 'weight', directed = False, missing_value = 0):
    if directed:
        u, v = edge
    else:
        u, v = edge
        if not G.has_edge(u, v):
            u, v = v, u  # 尝试反向
    data = G.get_edge_data(u, v, default={})
    return data.get(weight_key, missing_value) if data else missing_value

def aligned_edge_weights(G1, G2, weight_key='weight', directed=False):
    """
    返回两个图在边并集上的权重列表（顺序一致）
    """
    if directed:
        edges1 = set(G1.edges())
        edges2 = set(G2.edges())
    else:
        def norm(e):
            u, v = e
            return (u, v) if u <= v else (v, u)
        edges1 = {norm(e) for e in G1.edges()}
        edges2 = {norm(e) for e in G2.edges()}
    all_edges = sorted(edges1 | edges2)
    # def get_weight(G, edge):
    #    if directed:
    #        u, v = edge
    #    else:
    #        u, v = edge
    #        if not G.has_edge(u, v):
    #            u, v = v, u  # 尝试反向
    #    data = G.get_edge_data(u, v, default={})
    #    return data.get(weight_key, missing_value) if data else missing_value
    w1 = [get_weight(G1, e) for e in all_edges]
    w2 = [get_weight(G2, e) for e in all_edges]
    return all_edges, w1, w2

def OBMol_to_Graph(OB_mol):
    G = nx.Graph()
    # 添加带权重的边
    edges = []
    for bond in openbabel.OBMolBondIter(OB_mol):
        bond_idx = bond.GetIdx()
        begin_atom_idx = bond.GetBeginAtom().GetIdx()
        end_atom_idx = bond.GetEndAtom().GetIdx()
        bond_order = bond.GetBondOrder()
        edges.append( tuple([begin_atom_idx - 1, end_atom_idx - 1, bond_order]) )

    G.add_weighted_edges_from(edges)
    return G

def readmol(filename: str, fileformat: str):
    '''
    using openbabel to read the molecule
    input the file path and file format
    return a obmol object
    '''
    mol = openbabel.OBMol()
    obconversion = openbabel.OBConversion()
    obconversion.SetInFormat(fileformat)
    obconversion.ReadFile(mol,filename)
    return mol

def expand_subgraph_nx_elemental_ob(graph, initial_nodes, obmol, max_external_degree=1, batch=True, priority=None):
    """
    从初始节点集合开始，扩展子图直到所有边界节点的外部度数 ≤ max_external_degree。

    参数:
        graph: networkx.Graph，带权图（边需包含 'weight' 属性，可选）
        initial_nodes: iterable，初始子图的节点集合
        max_external_degree: int，允许的最大外部度数（默认2）
        batch: bool，是否批量添加（True：每次添加所有违规节点的所有外部邻居；
                           False：每次仅添加一个外部邻居，需配合priority使用）
        priority: str，当batch=False时，选择邻居的优先级
                  'weight'：优先添加权重最大的边对应的邻居
                  None：任意顺序（按邻居列表顺序）

    返回:
        set，扩展后的子图节点集合
    """
    # 使用集合方便快速成员判断
    subgraph = set(initial_nodes)

    # 辅助函数：计算当前子图中每个节点的外部度数
    def compute_external_degree():
        '''
        ext_deg = {}
        for node in subgraph:
            # 统计邻居中不在子图中的数量
            deg = sum(1 for nb in graph.neighbors(node) if nb not in subgraph)
            if deg > 0:
                ext_deg[node] = deg
        # print(ext_deg)
        return ext_deg
        '''
        ext_weight = {}
        for node in subgraph:
            # 累加所有外部邻居的边权重
            total = 0
            for nb in graph.neighbors(node):
                if nb not in subgraph:
                    # 只有C-C 键才做正常判断
                    if obmol.GetAtom(node + 1).GetAtomicNum() == 6 and obmol.GetAtom(nb + 1).GetAtomicNum() == 6:
                        w = graph[node][nb].get('weight', 1)   # 默认权重为1
                    else:
                        w = 2
                    # 非C-C键则通过加重权重，使得这个节点的外部连接度数超标，从而继续下一周期囊括
                    total += w
            if total > 0:
                ext_weight[node] = total
        return ext_weight

    # 批量添加模式
    if batch:
        while True:
            ext_deg = compute_external_degree()
            # 找出外部度数超限的节点
            violating = [n for n, d in ext_deg.items() if d > max_external_degree]
            if not violating:
                break
            # 收集所有需要添加的外部邻居
            to_add = set()
            for v in violating:
                for nb in graph.neighbors(v):
                    if nb not in subgraph:
                        to_add.add(nb)
            # 若无新节点可加（理论上不会发生），退出循环
            if not to_add:
                break
            subgraph.update(to_add)

    # 单步添加模式（每次只加一个邻居）
    else:
        while True:
            ext_deg = compute_external_degree()
            violating = [n for n, d in ext_deg.items() if d > max_external_degree]
            if not violating:
                break

            # 选择一个违规节点（这里简单取第一个，也可以优化为“最违规”的节点）
            v = violating[0]

            # 获取该节点的所有外部邻居及权重（如果有）
            external = [(nb, graph[v][nb].get('weight', 1)) for nb in graph.neighbors(v) if nb not in subgraph]
            if not external:
                continue  # 理论上不会发生

            # 根据优先级选择邻居
            if priority == 'weight':
                # 权重最大者
                nb, _ = max(external, key=lambda x: x[1])
            elif priority ==  'smallGraph':
                # 先侧基, 即断键后形成包含元素最少得连通图
                # 对于每个候选 nb，计算在移除当前子图后，nb 所在外部连通分量的大小
                # 注意：需要临时移除 subgraph 中的所有节点（但保留 nb 本身）
                # 使用 graph.subgraph 或手动构建剩余图
                # 为避免重复计算，提前构建剩余节点集合
                remaining_nodes = set(graph.nodes()) - subgraph
                # 在剩余节点诱导子图中，找到 nb 的连通分量
                # 注意：剩余节点可能包含多个连通分量，我们只需 nb 所在的那个
                # 使用 networkx 的 node_connected_component 方法，但需要指定子图
                # 更高效：构建一个只包含剩余节点的视图
                sub_remaining = graph.subgraph(remaining_nodes)
                # 计算每个候选 nb 所在分量的大小
                best_nb = None
                best_size = float('inf')
                for cand, _ in external:
                    # 如果 cand 不在 remaining_nodes 中？一定在，因为它是外部邻居
                    comp_size = len(nx.node_connected_component(sub_remaining, cand))
                    if comp_size < best_size:
                        best_size = comp_size
                        best_nb = cand
                    elif comp_size == best_size and best_nb is None:
                        best_nb = cand
                nb = best_nb

            else:
                # 默认取第一个
                nb = external[0][0]

            subgraph.add(nb)

    return subgraph

def expand_subgraph_nx_elemental_rd(graph, initial_nodes, rdmol , max_external_degree=1, batch=True, priority=None):
    """
    从初始节点集合开始，扩展子图直到所有边界节点的外部度数 ≤ max_external_degree。

    参数:
        graph: networkx.Graph，带权图（边需包含 'weight' 属性，可选）
        initial_nodes: iterable，初始子图的节点集合
        max_external_degree: int，允许的最大外部度数（默认2）
        batch: bool，是否批量添加（True：每次添加所有违规节点的所有外部邻居；
                           False：每次仅添加一个外部邻居，需配合priority使用）
        priority: str，当batch=False时，选择邻居的优先级
                  'weight'：优先添加权重最大的边对应的邻居
                  None：任意顺序（按邻居列表顺序）

    返回:
        set，扩展后的子图节点集合
    """
    # 使用集合方便快速成员判断
    subgraph = set(initial_nodes)
    
    # nodelist = []
    # for node in graph.nodes:
    #     nodelist.append(node)

    # print(subgraph)
    # print(nodelist)

    # 辅助函数：计算当前子图中每个节点的外部度数
    def compute_external_degree():
        '''
        ext_deg = {}
        for node in subgraph:
            # 统计邻居中不在子图中的数量
            deg = sum(1 for nb in graph.neighbors(node) if nb not in subgraph)
            if deg > 0:
                ext_deg[node] = deg
        # print(ext_deg)
        return ext_deg
        '''
        ext_weight = {}
        for node in subgraph:
            # 累加所有外部邻居的边权重
            total = 0
            for nb in graph.neighbors(node):
                if nb not in subgraph:
                    # 只有C-C 键才做正常判断
                    if rdmol.GetAtomWithIdx(node).GetAtomicNum() == 6 and rdmol.GetAtomWithIdx(nb + 1).GetAtomicNum() == 6:
                        w = graph[node][nb].get('weight', 1)   # 默认权重为1
                    else:
                        w = 2
                    # 非C-C键则通过加重权重，使得这个节点的外部连接度数超标，从而继续下一周期囊括
                    total += w
            if total > 0:
                ext_weight[node] = total
        return ext_weight

    # 批量添加模式
    if batch:
        while True:
            ext_deg = compute_external_degree()
            # 找出外部度数超限的节点
            violating = [n for n, d in ext_deg.items() if d > max_external_degree]
            if not violating:
                break
            # 收集所有需要添加的外部邻居
            to_add = set()
            for v in violating:
                for nb in graph.neighbors(v):
                    if nb not in subgraph:
                        to_add.add(nb)
            # 若无新节点可加（理论上不会发生），退出循环
            if not to_add:
                break
            subgraph.update(to_add)

    # 单步添加模式（每次只加一个邻居）
    else:
        while True:
            ext_deg = compute_external_degree()
            violating = [n for n, d in ext_deg.items() if d > max_external_degree]
            if not violating:
                break

            # 选择一个违规节点（这里简单取第一个，也可以优化为“最违规”的节点）
            v = violating[0]

            # 获取该节点的所有外部邻居及权重（如果有）
            external = [(nb, graph[v][nb].get('weight', 1)) for nb in graph.neighbors(v) if nb not in subgraph]
            if not external:
                continue  # 理论上不会发生

            # 根据优先级选择邻居
            if priority == 'weight':
                # 权重最大者
                nb, _ = max(external, key=lambda x: x[1])
            elif priority ==  'smallGraph':
                # 先侧基, 即断键后形成包含元素最少得连通图
                # 对于每个候选 nb，计算在移除当前子图后，nb 所在外部连通分量的大小
                # 注意：需要临时移除 subgraph 中的所有节点（但保留 nb 本身）
                # 使用 graph.subgraph 或手动构建剩余图
                # 为避免重复计算，提前构建剩余节点集合
                remaining_nodes = set(graph.nodes()) - subgraph
                # 在剩余节点诱导子图中，找到 nb 的连通分量
                # 注意：剩余节点可能包含多个连通分量，我们只需 nb 所在的那个
                # 使用 networkx 的 node_connected_component 方法，但需要指定子图
                # 更高效：构建一个只包含剩余节点的视图
                sub_remaining = graph.subgraph(remaining_nodes)
                # 计算每个候选 nb 所在分量的大小
                best_nb = None
                best_size = float('inf')
                for cand, _ in external:
                    # 如果 cand 不在 remaining_nodes 中？一定在，因为它是外部邻居
                    comp_size = len(nx.node_connected_component(sub_remaining, cand))
                    if comp_size < best_size:
                        best_size = comp_size
                        best_nb = cand
                    elif comp_size == best_size and best_nb is None:
                        best_nb = cand
                nb = best_nb

            else:
                # 默认取第一个
                nb = external[0][0]

            subgraph.add(nb)

    return subgraph


def run():
    args = parser_args()
    file_path = args.filename
    selection = args.selection
    sep = args.sep
    base = args.base
    output_pdb = args.output_pdb

    if sep is None:
        sep = " "
    if base is None:
        base = 0

    selected_atom_idx, charge = select_atoms(file_path, selection, output_pdb)
    selected_atom_idx = sorted(selected_atom_idx)
    selected_atom_idx_str = [ str(i + base) for i in selected_atom_idx ]
    if output_pdb is None:
        print(sep.join(selected_atom_idx_str))
        print("fragment charge calculated at multiplicity 1:", charge)


def output_RDmol(RDmol, atom_indices,  file_name):
    # atom_indices_tuple = tuple(atom_indices)
    bond_list = []
    for bond in RDmol.GetBonds():
        begin_atom_idx = bond.GetBeginAtom().GetIdx()
        end_atom_idx = bond.GetEndAtom().GetIdx()
        if begin_atom_idx in atom_indices and end_atom_idx in atom_indices:
            bond_list.append( bond.GetIdx() )

    try:
        selected_fragment = Chem.PathToSubmol(RDmol, bond_list)
    except Exception as e:
        print(f"extract fragment error, {e}")
        return False
    Chem.MolToPDBFile(selected_fragment, file_name)
    return True



def openmm_forcefield_struct_correction(file_path, selected_atom_idx): # this didn't help at all.
    forcefield = ForceField('amber14-all.xml', 'implicit/obc2.xml')
    openmm_mol = PDBFile(file_path)
    modeller = Modeller(openmm_mol.topology, openmm_mol.positions)
    # modeller.addHydrogens(forcefield)
    modeller.addHydrogens()

    restraint_force = CustomExternalForce("k * ( (x-x0)^2 + (y-y0)^2 + (z-z0)^2 )")
    restraint_force.addGlobalParameter("k", 1e7)
    restraint_force.addPerParticleParameter("x0")
    restraint_force.addPerParticleParameter("y0")
    restraint_force.addPerParticleParameter("z0")

    for atom, pos in zip( modeller.getTopology().atoms(), modeller.getPositions() ):
        if int(atom.index) in selected_atom_idx:
            restraint_force.addParticle( atom.index, [pos.x, pos.y, pos.z])


    system = forcefield.createSystem(modeller.topology, nonbondedMethod="PME")
    system.addForce(restraint_force)

    integrator = LangevinIntegrator(300*openmm_unit.kelvin, 1 / openmm_unit.picosecond, 0.002 * openmm_unit.picoseconds)
    platform = Platform.getPlatformByName('CPU')
    simulation = Simulation(modeller.topology, system, integrator, platform)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations = 50000)
    with open('tmp.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)


def select_atoms(file_path, selection, output_pdb):

    MDA_complex = mda.Universe(file_path)
    print("start selection", flush = True)
    selected_atoms = MDA_complex.select_atoms(selection)
    print("end selection", flush = True)
    selected_atom_idx = []
    
    for atom in selected_atoms:
        # print(atom.index, atom.element, atom.resname, atom.segid)
        selected_atom_idx.append(int(atom.index))



    print("RDKit reading", flush = True)
    RD_complex = Chem.MolFromPDBFile(file_path, removeHs = False, sanitize = True)
    if RD_complex is None:
        openmm_forcefield_struct_correction(file_path, selected_atom_idx)

    print("RDKit read", flush = True)
    OB_complex = None

    num_atoms = len(RD_complex.GetAtoms())
    obabel_max_num_threshold = 15000
    # one should avoid reading large pdb files using rdkit. Let's say when numatoms > 15000
    # TODO: write a rdkitmol to obabel mol converter
    if num_atoms < obabel_max_num_threshold:
        print("Openbabel reading", flush = True)
        OB_complex = readmol(file_path, "pdb")
        print("Openbabel read", flush = True)

    OB_G = None
    if OB_complex is not None:
        OB_G = OBMol_to_Graph(OB_complex)

    RD_G = RDMol_to_Graph(RD_complex)



    merged_G = None
    if OB_G is not None:
        merged_G = Merge_two_graphs_keep_high_BO(OB_G, RD_G)
    else:
        merged_G = RD_G

    expanded_atom_idx = selected_atom_idx
    if OB_complex is not None:
        expanded_atom_idx = expand_subgraph_nx_elemental_ob(merged_G, selected_atom_idx, OB_complex, batch=False, priority='smallGraph')
    else:
        expanded_atom_idx = expand_subgraph_nx_elemental_rd(merged_G, selected_atom_idx, RD_complex, batch=False, priority='smallGraph')

    # calc charge
    charge_list = calc_Gasteiger_Charge(RD_complex)
    # charge_list = calc_Formal_Charge(RD_complex)
    atomic_num_list = get_Atomic_Num(RD_complex)

    frac_charge_list = [ float(charge_list[i]) for i in expanded_atom_idx ]
    frac_partial_charge = sum(frac_charge_list)
    # print(frac_partial_charge)

    frac_atomic_num_list = [ int(atomic_num_list[i]) for i in expanded_atom_idx ]
    frac_atomic_num = sum(frac_atomic_num_list)


    frac_charge = judge_charge_multiplicity(frac_atomic_num, frac_partial_charge  )
    if output_pdb is not None:
        output_RDmol(RD_complex, expanded_atom_idx, output_pdb)

    return expanded_atom_idx, frac_charge



if __name__ == '__main__':
	run()





