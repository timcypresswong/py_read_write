# require openbabel via conda install conda-forge::openbabel
# require numpy and scipy
from openbabel import openbabel
import numpy as np
import itertools
import copy
from scipy.spatial.distance import pdist
import argparse


fullvalance = {
    1: 1,
    5: 3,
    6: 4,
    7: 3,
    8: 2,
    9: 1,
    15: 5,
    16: 2,
    17: 1,
    # 30: 4,
    35: 1,
    53: 1
    }


def parser_args():
	parser = argparse.ArgumentParser(description='This is the program to extract the pattern information from a file.')
	parser.add_argument('filename', help='input file name')
	parser.add_argument('outputname', help='output file name')
	args = parser.parse_args()
	return args

def writemol(mol, filename: str, fileformat: str):
    obconversion = openbabel.OBConversion()
    obconversion.SetOutFormat(fileformat)
    obconversion.WriteFile(mol, filename)

def readmol(filename: str, fileformat: str):
    mol = openbabel.OBMol()
    obconversion = openbabel.OBConversion()
    obconversion.SetInFormat(fileformat)
    obconversion.ReadFile(mol, filename)
    return mol

def run():
    args = parser_args()
    file_path = args.filename
    output = args.outputname

    mol = readmol(filename = file_path, fileformat = 'xyz')
    mol = addH(mol)
    writemol(mol, filename = output, fileformat = 'xyz')


def are_atoms_collinear_vectors(coords, threshold=0.1):
    """
    使用向量叉积法判断四个原子是否共线

    参数:
    coords: 4个原子的坐标，形状为(4, 3)的数组
    threshold: 最大允许偏离距离（单位：Å）

    返回:
    bool: 是否共线
    float: 最大偏离距离
    """
    points = np.array(coords)

    if points.shape != (4, 3) and points.shape != (3, 3) :
        raise ValueError("需要4个原子或3个原子的三维坐标")

    p1 = None
    p2 = None
    p3 = None
    p4 = None

    if points.shape == (4, 3):
    # 选择前三个点作为参考
        p1, p2, p3, p4 = points
    else:
        p1, p2, p3 = points
    # 计算两个向量
    v1 = p2 - p1
    v2 = p3 - p1

    # 计算叉积的模长（平行四边形的面积）
    cross_product = np.cross(v1, v2)
    cross_norm = np.linalg.norm(cross_product)

    # 如果叉积模长为0或很小，则三点共线
    if cross_norm < threshold:
        # 三点为输入时, 直接返回true
        if p4 is None:
            return True, cross_norm

        # 三点共线，检查第四点
        # 计算直线方向向量
        line_direction = v1 / np.linalg.norm(v1) if np.linalg.norm(v1) > 0 else v2 / np.linalg.norm(v2)

        # 计算第四点到直线的距离
        # 点到直线的距离公式: d = |(p4-p1) × direction| / |direction|
        # 由于direction是单位向量，|direction|=1
        v4 = p4 - p1
        distance = np.linalg.norm(np.cross(v4, line_direction))

        is_collinear = distance < threshold
        return is_collinear, distance
    else:
        # 前三点不共线，四个点不可能共线
        return False, cross_norm

def are_atoms_coplanar(coords, threshold=0.1):
    """
    判断原子是否近似共面

    参数:
    coords: 6个原子的坐标，形状为(6, 3)的数组
    threshold: 最大允许偏离距离（单位：Å），默认0.1Å

    返回:
    bool: 是否共面
    float: 最大偏离距离
    """
    # 确保输入是numpy数组
    points = np.array(coords)

    # if points.shape != (6, 3):
    if points.shape[0] < 4:

        raise ValueError("需要至少4个原子的三维坐标")

    # 1. 中心化数据
    centroid = np.mean(points, axis=0)
    centered_points = points - centroid

    # 2. 使用SVD找到最佳拟合平面
    # 协方差矩阵
    cov_matrix = np.dot(centered_points.T, centered_points) / (len(points) - 1)

    # 奇异值分解
    U, S, Vt = np.linalg.svd(cov_matrix)

    # 最小奇异值对应的特征向量即为平面的法向量
    normal_vector = Vt[2, :]

    # 3. 计算每个点到平面的距离
    distances = np.abs(np.dot(centered_points, normal_vector))

    # 4. 计算平面方程
    # 平面方程: ax + by + cz + d = 0
    a, b, c = normal_vector
    d = -np.dot(normal_vector, centroid)

    # 5. 判断是否共面
    max_distance = np.max(distances)
    is_coplanar = max_distance < threshold

    return is_coplanar, max_distance, (a, b, c, d)

def inRing_AtomIdx(mol, ring):
    atomidxlist = []
    #print("ring path:", ring._path)
    for obatom in openbabel.OBMolAtomIter(mol):
        if ring.IsInRing(obatom.GetIdx()):
            atomidxlist.append( obatom.GetIdx() )
    return atomidxlist

def is_aromatic_positional(mol, atomlist) -> bool:
    augmented_atom_list = copy.deepcopy(atomlist)
    for idx in atomlist:
        for neighatom in openbabel.OBAtomAtomIter(mol.GetAtom(idx)):
            augmented_atom_list.append(neighatom.GetIdx())
    augmented_atom_list = list(set(augmented_atom_list))
    atom_coord = []
    for idx in augmented_atom_list:
        atom = mol.GetAtom(idx)
        atom_coord.append( [atom.GetX(), atom.GetY(), atom.GetZ() ]  )
    atom_coord = np.array(atom_coord)
    if atom_coord.shape[0] > 3:
        if are_atoms_coplanar(atom_coord, threshold=0.18)[0]:
            return True
    # print("the ring is not coplanar enough:", atomlist, are_atoms_coplanar(atom_coord, threshold=0.15)[1])
    return False

def is_aromatic_electron(mol, atomlist):
    count_electron = 0
    N_idx = []                              # the index of N that can be work around
    num_double_bond = 0
    non_bond_N_num = 0
    for idx in atomlist:
        Atom = mol.GetAtom(idx)
        AtomicNum = Atom.GetAtomicNum()
        if AtomicNum == 6:
            if Atom.GetExplicitValence() < 4:
                count_electron += 1
                num_double_bond += 1
        if AtomicNum == 8:
            count_electron += 2
        if AtomicNum == 16:
            count_electron += 2
        if AtomicNum == 7:
            if Atom.GetExplicitValence() == 3:
                count_electron += 2
            if Atom.GetExplicitValence() == 2:
                count_electron += 1
                N_idx.append(idx)
                num_double_bond += 1

    if count_electron > 15:
        # dense ring (possibly more than 15 like 芘 reveals an anti Huckel rule, as long as they are planar enough, we consider them aromatic)
        num_double_bond = int(num_double_bond / 2)
        return  True, N_idx, num_double_bond, non_bond_N_num
    elif (count_electron % 4 - 2) == 0:
        num_double_bond = int(num_double_bond / 2)
        return  True, N_idx, num_double_bond, non_bond_N_num
    elif len(N_idx) >= np.abs(count_electron % 4 - 2 ):
        num_double_bond = num_double_bond - np.abs(count_electron % 4 - 2 )
        num_double_bond = int(num_double_bond / 2)
        non_bond_N_num =   np.abs(count_electron % 4 - 2 )
        return True, N_idx, num_double_bond, non_bond_N_num

    # print("electron in this ring is not considered as aromatic", count_electron)
    return False, N_idx, 0, non_bond_N_num

def max_distance_sum_bruteforce(coords, M):
    """
    暴力枚举所有组合，选择距离和最大的M个点

    参数:
    coords: (N, 3) numpy数组，氮原子坐标
    M: 要选择的原子数

    返回:
    best_indices: 选择的M个原子的索引
    max_sum: 最大距离和
    """
    N = len(coords)

    if M > N:
        raise ValueError(f"M({M})不能大于N({N})")

    if M == 0:
        # 如果只选1个或0个点，距离和为0
        return [], 0.0

    if M == 1:
        # 如果只选1个或0个点，距离和为0
        return [0], 0.0


    best_indices = None
    max_sum = -float('inf')

    # 遍历所有可能的组合
    for indices in itertools.combinations(range(N), M):
        # 获取选中的坐标
        selected_coords = coords[list(indices)]

        # 计算所有点对之间的距离
        # 方法1：使用pdist计算所有点对距离
        distances = pdist(selected_coords)

        # 计算距离总和
        dist_sum = np.sum(distances)

        # 更新最大值
        if dist_sum > max_sum:
            max_sum = dist_sum
            best_indices = indices

    return list(best_indices), max_sum

def find_double_edges(edges, M):
    """
    在图中选择M条边升级为双键，使所有节点度数≤3

    参数:
        edges: 边列表 [[u1,v1], [u2,v2], ...]
        M: 要升级的边数

    返回:
        选择的边索引列表，或None表示无解
    """
    N = len(edges)

    # # 初始度数计算
    # def compute_degrees(selected_edges_indices):
    #     degree = {}
    #     # 初始：所有边都是单键
    #     for u, v in edges:
    #         degree[u] = degree.get(u, 0) + 1
    #         degree[v] = degree.get(v, 0) + 1

    #     # 为选中的边升级（增加度数）
    #     for idx in selected_edges_indices:
    #         u, v = edges[idx]
    #         degree[u] += 1  # 单键变双键，度数+1
    #         degree[v] += 1

    #     return degree

    # # 检查是否满足度数约束
    # def is_valid(selected_edges_indices):
    #     degree = compute_degrees(selected_edges_indices)
    #     # print(selected_edges_indices, degree)
    #     return all(d <= 3 for d in degree.values())
        # 计算度数变化量
  # 计算初始度数（所有边都是单键）
    def compute_initial_degrees():
        initial_degrees = {}
        for u, v in edges:
            initial_degrees[u] = initial_degrees.get(u, 0) + 1
            initial_degrees[v] = initial_degrees.get(v, 0) + 1
        return initial_degrees

    # 计算与初始状态相比的变化量
    def compute_change_amounts(selected_indices, initial_degrees):
        """计算每个节点相对于初始状态的变化量"""
        # 计算升级后的度数
        upgraded_degrees = initial_degrees.copy()

        # 对于选中的边，度数增加1（因为单键变双键）
        for idx in selected_indices:
            u, v = edges[idx]
            upgraded_degrees[u] += 1
            upgraded_degrees[v] += 1

        # 计算变化量 = 升级后度数 - 初始度数
        change_amounts = {}
        for node, upgraded_deg in upgraded_degrees.items():
            initial_deg = initial_degrees.get(node, 0)
            change_amount = upgraded_deg - initial_deg
            if change_amount != 0:  # 只记录有变化的节点
                change_amounts[node] = change_amount

        return change_amounts

    # 计算初始度数（所有边都是单键时的度数）
    initial_degrees = compute_initial_degrees()

    # 检查是否满足变化量约束
    def is_valid(selected_indices):
        change_amounts = compute_change_amounts(selected_indices, initial_degrees)
        # print(selected_indices, change_amounts)
        # 关键：变化量必须 < 2（不允许大于等于2）
        return all(change < 2 for change in change_amounts.values())
    # 回溯搜索
    result = []

    def backtrack(start, selected):
        nonlocal result

        # 如果已经找到解，直接返回
        if result:
            return

        # 如果已经选择了M条边，检查是否有效
        if len(selected) == M:
            if is_valid(selected):
                result = selected.copy()
            return

        # 如果剩余边不够，剪枝
        if len(selected) + (N - start) < M:
            return

        # 尝试选择当前边
        selected.append(start)
        if is_valid(selected):  # 只有有效才继续
            backtrack(start + 1, selected)
        selected.pop()

        # 尝试不选择当前边
        backtrack(start + 1, selected)

    backtrack(0, [])
    return result if result else None

def generate_ring_information_for_aromatic_calculation(mol, atomList):
    aromatic_list = []
    aromatic_Nlist = []
    num_double_bond_list = []
    non_bond_N_num_list = []
    for ringlist in atomList:
        if is_aromatic_positional(mol, ringlist):
            aromatic_elec, Nlist, num_double_bond, non_bond_N_num = is_aromatic_electron(mol, ringlist)
            aromatic_list.append(aromatic_elec)
            aromatic_Nlist.append(Nlist)
            num_double_bond_list.append(num_double_bond)
            non_bond_N_num_list.append(non_bond_N_num)
        else:
            aromatic_list.append(False)
            aromatic_Nlist.append([])
            num_double_bond_list.append(0)
            non_bond_N_num_list.append(0)


    provide_2_electron_N_idx_list = [ [] for i in aromatic_list ]
    for ring_idx in range(len(aromatic_list)):
        if aromatic_list[ring_idx]:
            N_coord = []
            to_remove_N_idx = []
            for idx in aromatic_Nlist[ring_idx]:
                atom = mol.GetAtom(idx)
                N_coord.append( [atom.GetX(), atom.GetY(), atom.GetZ() ]  )

            N_coord = np.array(N_coord)

            to_remove_N_idx = max_distance_sum_bruteforce(N_coord, non_bond_N_num_list[ring_idx])[0]
            # print(N_coord,non_bond_N_num, to_remove_N_idx)
            provide_2_electron_N_idx = [idx for i, idx in enumerate(aromatic_Nlist[ring_idx]) if i in to_remove_N_idx ]
            provide_2_electron_N_idx_list[ring_idx] = provide_2_electron_N_idx
    return aromatic_list, num_double_bond_list, non_bond_N_num_list, provide_2_electron_N_idx_list

def merge_overlapping_sets(sets):
    '''
    合并列表中有重叠元素的集合。
    参数 sets: 列表，每个元素是可转换为集合的类型（如元组、列表、集合）。
    返回: 合并后的集合列表。
    '''
    # 将所有元素转换为集合
    sets = [set(s) for s in sets]
    n = len(sets)

    # 并查集初始化
    parent = list(range(n))

    def find(x):
        # 路径压缩
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[ry] = rx

    # 构建元素 -> 索引列表的映射
    elem_to_indices = {}
    for idx, s in enumerate(sets):
        for elem in s:
            elem_to_indices.setdefault(elem, []).append(idx)

    # 合并包含相同元素的集合的索引
    for indices in elem_to_indices.values():
        if len(indices) > 1:
            first = indices[0]
            for other in indices[1:]:
                union(first, other)

    # 将同一根节点的所有集合合并
    groups = {}
    for idx, s in enumerate(sets):
        root = find(idx)
        if root not in groups:
            groups[root] = set()
        groups[root].update(s)

    return list(groups.values())

def correct_bond_order(mol):

    # purge
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if bond.IsInRing():


            if obatom1.GetAtomicNum() == 6 and  (obatom2.GetAtomicNum() == 6 or obatom2.GetAtomicNum() == 7 or obatom2.GetAtomicNum() == 8) :
                bond.SetBondOrder(1)
            if obatom2.GetAtomicNum() == 6 and  (obatom1.GetAtomicNum() == 6 or obatom1.GetAtomicNum() == 7 or obatom1.GetAtomicNum() == 8) :
                bond.SetBondOrder(1)


    # C - O bond
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 8 and obatom1.GetExplicitValence() == 1 and obatom2.GetExplicitValence() < fullvalance[ obatom2.GetAtomicNum()]:
            atom_coords = []

            for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
            for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
            atom_coords = np.array(atom_coords)
            if bond.GetLength() < 1.25:
                if atom_coords.shape == (4, 3):
                    if are_atoms_coplanar(atom_coords)[0]:
                        bond.SetBondOrder(2)

                if atom_coords.shape == (6, 3):
                    if are_atoms_coplanar(atom_coords)[0]:
                        bond.SetBondOrder(2)
                if atom_coords.shape == (3, 3):
                    bond.SetBondOrder(2)

        elif obatom2.GetAtomicNum() == 8 and obatom2.GetExplicitValence() == 1 and obatom1.GetExplicitValence() < fullvalance[ obatom1.GetAtomicNum()]:
            atom_coords = []

            for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
            for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
            atom_coords = np.array(atom_coords)
            if bond.GetLength() < 1.25:

                if atom_coords.shape == (4, 3):
                    if are_atoms_coplanar(atom_coords)[0]:
                        bond.SetBondOrder(2)

                if atom_coords.shape == (6, 3):
                    if are_atoms_coplanar(atom_coords)[0]:
                        bond.SetBondOrder(2)
                if atom_coords.shape == (3, 3):
                    bond.SetBondOrder(2)

    CCdoublebond_added = []

    SmallestRingList = []
    for obring in openbabel.OBMolRingIter(mol):
        atomlist =  inRing_AtomIdx(mol, obring)
        SmallestRingList.append(atomlist)
        # print(obring.Size(), atomlist)

    # molGraph = mol_to_graph(mol)
    # SelfList = find_all_chordless_cycles(molGraph)
    # print("network x implement:", SelfList)

    aromatic_list, num_double_bond_list, non_bond_N_num_list, provide_2_electron_N_idx_list = generate_ring_information_for_aromatic_calculation(mol, SmallestRingList)


    # combine the no chord rings to a larger one if the rings are connected to each other
    CombinedRingList = []

    flag_append_list = [True for x in range(len(SmallestRingList))]
    ring_connection_list = []
    for i, ringlist in enumerate(SmallestRingList):
        for j in range(i + 1, len(SmallestRingList)):
            intersection = set.intersection( set(SmallestRingList[i]), set(SmallestRingList[j])  )
            # union = set.union( set(SmallestRingList[i]), set(SmallestRingList[j])  )
            if len(intersection) == 2 and aromatic_list[i] == True and  aromatic_list[j] == True:
                ring_connection_list.append( set([i,j])   )
                # CombinedRingList.append(list(union))
                flag_append_list[i] = False
                flag_append_list[j] = False
        if flag_append_list[i]:
            CombinedRingList.append(SmallestRingList[i])
    merged_connection_list = merge_overlapping_sets(ring_connection_list)
    for ringidxlist in merged_connection_list:
        atom_in_bigring = []
        for ringidx in ringidxlist:
            atom_in_bigring = set.union( set(atom_in_bigring), set( SmallestRingList[ringidx] )  )
        CombinedRingList.append( list(atom_in_bigring)  )

    #print(CombinedRingList)
    aromatic_list, num_double_bond_list, non_bond_N_num_list, provide_2_electron_N_idx_list = generate_ring_information_for_aromatic_calculation(mol, CombinedRingList)

    #print(aromatic_list)
    #print(num_double_bond_list)
    #print(non_bond_N_num_list)
    #print(provide_2_electron_N_idx_list)

    for ringidx in range(len(aromatic_list)):
        if aromatic_list[ringidx] == True:
            possible_modified_bond_list = []
            for bond in openbabel.OBMolBondIter(mol):
                obatom1 = bond.GetBeginAtom()
                obatom2 = bond.GetEndAtom()
                if obatom1.GetIdx() in CombinedRingList[ringidx] and obatom2.GetIdx() in CombinedRingList[ringidx]:
                    if (obatom1.GetExplicitValence() < fullvalance[obatom1.GetAtomicNum()]) and ( obatom2.GetExplicitValence() < fullvalance[obatom2.GetAtomicNum()] ):
                        if obatom1.GetIdx() not in provide_2_electron_N_idx_list[ringidx] and obatom2.GetIdx() not in provide_2_electron_N_idx_list[ringidx]:
                            possible_modified_bond_list.append( [obatom1.GetIdx() , obatom2.GetIdx() ] )
            # print("atom idx in CombinedRingList:", CombinedRingList[ringidx])
            # print("possible modified bond list:", possible_modified_bond_list)
            # print("number of double bond to add:", num_double_bond_list[ringidx])
            double_edges_idx = find_double_edges(possible_modified_bond_list,num_double_bond_list[ringidx] )
            # print("double edges idx:", double_edges_idx)
            double_edges = [edges for i, edges in enumerate(possible_modified_bond_list) if i in double_edges_idx]
            # print("double edges:", double_edges)
            for pair in double_edges:
                mol.GetBond(pair[0], pair[1]).SetBondOrder(2)


    # SP2-SP2 C=C-Heter
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 6 and obatom2.GetAtomicNum() == 6:
            atom_coords = []
            atom_Atomic_nums = []
            if obatom1.GetExplicitValence() == 3 and obatom2.GetExplicitValence() == 3:
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )

                    atom_Atomic_nums.append( neighbouratom.GetAtomicNum() )
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                    atom_Atomic_nums.append( neighbouratom.GetAtomicNum() )

                atom_coords = np.array(atom_coords)
                if atom_coords.shape == (6, 3) and atom_Atomic_nums != [6,6,6,6,6,6]:
                    if are_atoms_coplanar(atom_coords,threshold= 0.15)[0]:
                        bond.SetBondOrder(2)
                        CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )

    # SP2-SP2 C=C-C
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 6 and obatom2.GetAtomicNum() == 6:
            atom_coords = []
            # atom_coords.append( [obatom1.GetX(), obatom1.GetY(), obatom1.GetZ() ] )
            # atom_coords.append( [obatom2.GetX(), obatom2.GetY(), obatom2.GetZ() ] )
            if obatom1.GetExplicitValence() == 3 and obatom2.GetExplicitValence() == 3:
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                # atom_coords = list(set(atom_coords))
                atom_coords = np.array(atom_coords)
                # print(atom_coords)
                if atom_coords.shape == (6, 3):
                    if are_atoms_coplanar(atom_coords,threshold= 0.15)[0]:
                        bond.SetBondOrder(2)
                        CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )


    # triple bond
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 6 or obatom2.GetAtomicNum() == 6:
            atom_coords = []
            if obatom1.GetExplicitValence() == 2 or obatom2.GetExplicitValence() == 2:
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                atom_coords = np.array(atom_coords)
                if atom_coords.shape == (4, 3) or atom_coords.shape == (3, 3):
                    if are_atoms_collinear_vectors(atom_coords)[0]:
                        bond.SetBondOrder(3)

    # C-N
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 7 and obatom2.GetAtomicNum() == 6:
            if (obatom1.GetExplicitValence() < fullvalance[obatom1.GetAtomicNum()]) and (obatom2.GetExplicitValence() < fullvalance[obatom2.GetAtomicNum()] - 1) and bond.GetLength() < 1.39:
                bond.SetBondOrder(2)
                CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )
        if obatom2.GetAtomicNum() == 7 and obatom1.GetAtomicNum() == 6:
            if (obatom2.GetExplicitValence() < fullvalance[obatom2.GetAtomicNum()]) and (obatom1.GetExplicitValence() < fullvalance[obatom1.GetAtomicNum()] - 1) and bond.GetLength() < 1.39:
                bond.SetBondOrder(2)
                CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )

    # SP2-SP2 C=C-C
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 6 and obatom2.GetAtomicNum() == 6:
            atom_coords = []
            if obatom1.GetExplicitValence() == 3 or obatom2.GetExplicitValence() == 3:
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                atom_coords = np.array(atom_coords)
                if atom_coords.shape == (5, 3):
                    if are_atoms_coplanar(atom_coords, threshold= 0.15)[0]:
                        if bond.IsInRing():
                            continue
                            # bond.SetBondOrder(2)
                            # CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )
                        elif bond.GetLength() < 1.35:
                            bond.SetBondOrder(2)
                            CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )

    # SP2-SP2 C=C-C
    for bond in openbabel.OBMolBondIter(mol):
        obatom1 = bond.GetBeginAtom()
        obatom2 = bond.GetEndAtom()
        if obatom1.GetAtomicNum() == 6 and obatom2.GetAtomicNum() == 6:
            atom_coords = []
            if obatom1.GetExplicitValence() == 2 or obatom2.GetExplicitValence() == 2:
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom1):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                for neighbouratom  in openbabel.OBAtomAtomIter(obatom2):
                    atom_coords.append( [neighbouratom.GetX(), neighbouratom.GetY(), neighbouratom.GetZ() ] )
                atom_coords = np.array(atom_coords)
                if atom_coords.shape == (4, 3):
                    if are_atoms_coplanar(atom_coords, threshold= 0.15)[0]:
                        if bond.IsInRing():
                            continue
                            # bond.SetBondOrder(2)
                            # CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )
                        elif bond.GetLength() < 1.35:
                            bond.SetBondOrder(2)
                            CCdoublebond_added.append(  [obatom1.GetIdx(), obatom2.GetIdx() ]  )

    C_index_unique = []
    for pair in CCdoublebond_added:
        C_index_unique.append( pair[0] )
        C_index_unique.append( pair[1] )
    C_index_unique = list(set(C_index_unique))
    # C_count_dict = { x: 0 for x in C_index_unique }
    C_count_dict = { x: 0 for x in range(1, mol.NumAtoms() + 1) }

    for pair in CCdoublebond_added:
        C_count_dict[pair[0]] += 1
        C_count_dict[pair[1]] += 1
    # print(CCdoublebond_added)

    for key, value in C_count_dict.items():
        if value > 1:
            for pair in CCdoublebond_added :
                if key in pair and mol.GetAtom(pair[0]).GetAtomicNum() == 6 and mol.GetAtom(pair[1]).GetAtomicNum() == 6 : #
                    other = pair[1]
                    if pair[1] == key:
                        other = pair[0]
                    atom = mol.GetAtom( other )

                    oldbond = mol.GetBond( pair[0], pair[1] )
                    oldbond.SetBondOrder(1)
                    C_count_dict[key] -= 1
                    C_count_dict[other] -= 1
                    # print("setting", key, other)
                    for neighbouratom  in openbabel.OBAtomAtomIter(atom):
                        neighIdx = neighbouratom.GetIdx()
                        neighAtomicNum = neighbouratom.GetAtomicNum()
                        # print(key, other, neighIdx)
                        if (neighIdx not in pair) and (neighAtomicNum == 6) and (neighbouratom.GetExplicitValence() <= neighbouratom.GetExplicitDegree()) :
                            # print(key, other, neighIdx, "yes")
                            newbond = mol.GetBond( other , neighIdx )
                            if newbond.IsInRing():
                                continue
                                newbond.SetBondOrder(2)
                                C_count_dict[other] += 1
                                C_count_dict[neighIdx] += 1

                                # try:
                                #     C_count_dict[neighIdx] += 1
                                # except KeyError:
                                #     C_count_dict[neighIdx] = 1

                            elif newbond.GetLength() < 1.35:
                                newbond.SetBondOrder(2)
                                C_count_dict[other] += 1
                                C_count_dict[neighIdx] += 1

                                # try:
                                #     C_count_dict[neighIdx] += 1
                                # except KeyError:
                                #     C_count_dict[neighIdx] = 1

                    break

    return mol

def addH(mol):
    mol = correct_bond_order(mol)

    for obatom in openbabel.OBMolAtomIter(mol):
        valancediff = fullvalance[ obatom.GetAtomicNum() ] - obatom.GetTotalValence()
        if (obatom.GetImplicitHCount() + valancediff > 0):
            obatom.SetImplicitHCount( obatom.GetImplicitHCount() + valancediff )
    mol.AddHydrogens()
    return mol


if __name__ == '__main__':
	run()
