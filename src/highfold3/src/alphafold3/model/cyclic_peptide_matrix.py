import numpy as np
from numpy import ndarray


#def calc_cyclic_matrix_af3(residue_index,
#                           head_to_tail: bool = False,
#                           list_cid_ss: list[list[int]] = [],
#                           token_wise=True,
#                           signal_bug_fix=True, ):
#  # m = np.arange(len(residue_index))
#  # offset = m[:, None] - m[None, :]
#  # mn = np.stack([m, m + len(residue_index)], -1)
#  # c_offset = np.abs(mn[:, None, :, None] - mn[None, :, None, :]).min((2, 3))
#  # if signal_bug_fix:
#  #   a = c_offset < np.abs(offset)
#  #   c_offset[a] = -c_offset[a]
#  # final_offset = c_offset * np.sign(offset)
#  # return final_offset  
#  if not head_to_tail and (list_cid_ss is None or list_cid_ss == [] or list_cid_ss == [[]]):
#    return [[0]]
#  index_residue_1 = np.where(residue_index == 1)[0]
#  index_cid = 0
#  for i in range(len(list_cid_ss)):
#    if list_cid_ss[i][0] > index_cid:
#      index_cid = list_cid_ss[i][0]
#  if index_cid == 0:
#    return [[0]]
#  if len(index_residue_1) < index_cid:
#    print('error input chain id in disulfide bridges information')
#    return [[0]]
#  m = np.arange(len(residue_index))
#  if index_cid < len(index_residue_1):
#    m = np.arange(index_residue_1[index_cid])
#  if not token_wise:
#    m = residue_index[:len(m)]
#  offset = m[:, None] - m[None, :]
#  for i in range(len(list_cid_ss)):
#    index_cid_ss = list_cid_ss[i]
#    res_index_0 = index_residue_1[index_cid_ss[0] - 1]
#    n_res = index_residue_1[index_cid_ss[0]] if len(index_residue_1) > index_cid_ss[0] else len(residue_index)
#    n_res -= res_index_0
#    c1 = [ind - 1 for ind in index_cid_ss[1::2]]
#    c2 = [ind - 1 for ind in index_cid_ss[2::2]]
#    matrix = calc_cyclic_matrix_monomer(n_res, c1, c2, head_to_tail=head_to_tail)
#    offset[res_index_0:res_index_0 + n_res, res_index_0:res_index_0 + n_res] = matrix
#  # print(offset)
#  return offset  

def calc_cyclic_matrix_af3(residue_index,
                           head_to_tail: bool = False,
                           list_cid_ss: list[list[int]] = [],
                           token_wise=True,
                           signal_bug_fix=True):
  if not head_to_tail and (list_cid_ss is None or list_cid_ss == [] or list_cid_ss == [[]]):
    return [[0]]
  if len(list_cid_ss) > 1:
    return [[0]]
  index_residue_1 = np.where(residue_index == 1)[0]
  n_res = index_residue_1[1] if len(index_residue_1) > 1 else len(residue_index)
  # m = np.arange(len(n_res))
  # offset = m[:, None] - m[None, :]
  c1 = []
  c2 = []
  if len(list_cid_ss) > 0:
    c1 = [ind - 1 for ind in list_cid_ss[0][1::2]]
    c2 = [ind - 1 for ind in list_cid_ss[0][2::2]]
  offset = calc_cyclic_matrix_monomer(n_res, c1, c2, head_to_tail=head_to_tail)
  return offset

def calc_cyclic_matrix_monomer(n_aa, c1, c2, head_to_tail=False):
  """
  :param n_aa: number of amino acid residues
  :param c1: list of Cys residues
  :param c2: list of the corresponding Cys residues to c1, the same length with c1
  :param head_to_tail:
  :return:
  """
  if len(c1) != len(c2):
    return []  
  # init adjacency matrix
  matrix = np.zeros((n_aa, n_aa)) + n_aa
  for i in range(n_aa):
    matrix[i][i] = 0  
  # linear peptide connection
  for i in range(n_aa - 1):
    matrix[i][i + 1] = 1
    matrix[i + 1][i] = 1  
  # nc connection
  if head_to_tail:
    matrix[0][n_aa - 1] = 1
    matrix[n_aa - 1][0] = 1  
  # ss connection
  for i in range(len(c1)):
    matrix[c1[i]][c2[i]] = 1
    matrix[c2[i]][c1[i]] = 1  
  # get the shortest path
  matrix = calc_shortest_path(matrix)  
  for i in range(matrix.shape[0]):
    for j in range(i + 1, matrix.shape[0], 1):
      matrix[i][j] *= -1  
  return matrix


def calc_shortest_path(matrix):
  """
  Floyd algorithm to find the shortest path
  :param matrix:
  :return:
  """
  path = np.zeros_like(matrix)
  for i in range(matrix.shape[0]):
    path[i] = [j for j in range(matrix.shape[0])]
  # print(path)
  # print()
  # print(matrix)  
  for m in range(matrix.shape[0]):
    for i in range(matrix.shape[0]):
      for j in range(matrix.shape[0]):
        if matrix[i][m] + matrix[m][j] < matrix[i][j]:
          matrix[i][j] = matrix[i][m] + matrix[m][j]
          path[i][j] = m
  # print()
  # print(path)
  # print()
  # print(matrix)  
  return matrix


def construct_cyclic_offset_matrix_v20250306(n_res, signal_bug_fix=True,):
  i = np.arange(n_res)
  offset = i[:, None] - i[None, :]
  ij = np.stack([i, i + n_res], -1)
  c_offset = np.abs(ij[:, None, :, None] - ij[None, :, None, :]).min((2, 3))
  if signal_bug_fix:
    a = c_offset < np.abs(offset)
    c_offset[a] = -c_offset[a]
  final_offset = c_offset * np.sign(offset)
  return final_offset



def calc_cyclic_token_matrices(
    residue_index, token_index, asym_id, chain_type, res_name, ccd,
    *, protein_chain_type, head_to_tail=False, list_cid_ss=None,
):
  """Map the upstream residue graph onto AF3 tokens in the first peptide chain.

  A modified residue can occupy several tokens. Chain boundaries come from
  asym_id, while disulfide positions remain one-based residue positions.
  Within-residue atom offsets retain the historical AF3 token encoding.
  """
  bridges = [] if list_cid_ss is None else list_cid_ss
  if not isinstance(bridges, list):
    raise ValueError('disulfide_chain_res must be a JSON list of lists')
  if not head_to_tail and not bridges:
    zero = np.zeros((1, 1), dtype=np.int32)
    return zero, zero.copy()
  if len(residue_index) == 0 or chain_type[0] != protein_chain_type:
    raise ValueError('The first chain must be a peptide for HighFold3 cyclization')
  positions = np.flatnonzero(asym_id == asym_id[0])
  if not np.array_equal(positions, np.arange(len(positions))):
    raise ValueError('The first peptide chain must occupy the initial token block')
  first_residues = np.asarray(residue_index)[positions]
  residues = np.unique(first_residues)
  if not np.array_equal(residues, np.arange(1, len(residues) + 1)):
    raise ValueError('Peptide residue indices must be consecutive and start at 1')
  if len(residues) < 2:
    raise ValueError('Cyclization requires at least two peptide residues')
  if len(bridges) > 1:
    raise ValueError('This HighFold3 deployment supports one cyclic peptide chain')
  left, right = [], []
  if bridges:
    entry = bridges[0]
    if not isinstance(entry, list) or any(type(x) is not int for x in entry):
      raise ValueError('Disulfide entries must contain integer chain and residue indices')
    if len(entry) < 3 or len(entry) % 2 != 1 or entry[0] != 1:
      raise ValueError('Use [[1, residue1, residue2, ...]] for the first peptide chain')
    used = set()
    for a, b in zip(entry[1::2], entry[2::2]):
      if not 1 <= a <= len(residues) or not 1 <= b <= len(residues) or a == b:
        raise ValueError(f'Invalid disulfide residue pair: {a}, {b}')
      if a in used or b in used:
        raise ValueError('A cysteine cannot occur in more than one disulfide bridge')
      for residue in (a, b):
        name = res_name[positions[np.flatnonzero(first_residues == residue)[0]]]
        component = ccd.get(str(name))
        if component is None or 'SG' not in component.get('_chem_comp_atom.atom_id', ()):
          raise ValueError(f'Disulfide residue {residue} ({name}) has no SG atom')
      used.update((a, b))
      left.append(a - 1)
      right.append(b - 1)
  graph = np.asarray(
      calc_cyclic_matrix_monomer(len(residues), left, right, head_to_tail),
      dtype=np.int32,
  )
  residue_matrix = graph[(first_residues - 1)[:, None], (first_residues - 1)[None, :]]
  first_tokens = np.asarray(token_index)[positions]
  token_matrix = first_tokens[:, None] - first_tokens[None, :]
  return residue_matrix, np.asarray(token_matrix, dtype=np.int32)
