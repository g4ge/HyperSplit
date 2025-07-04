import time
import config
import argparse
import pandas as pd
from timeout_decorator import timeout, TimeoutError


def parse_args():
    parser = argparse.ArgumentParser()
    
    parser.add_argument("--x", type=str, help="CSV filename for table x")
    parser.add_argument("--y", type=str, help="CSV filename for table y")
    
    parser.add_argument("--xh", type=int, default=0, help="Number of header rows to skip in table x")
    parser.add_argument("--yh", type=int, default=0, help="Number of header rows to skip in table y")
    
    parser.add_argument("--d", type=float, default=None, help="Delta value for tolerance-based pruning")
    parser.add_argument("--n", type=int, default=None, help="Number of search tree levels to apply tolerance-based pruning")
    
    return parser.parse_args()


class HpgNode:
    def __init__(self, value, row, col, deg):
        self.value = value
        self.row = row
        self.col = col
        self.deg = deg


class NodePair:
    def __init__(self, v, w):
        self.v = v
        self.w = w


class Bidomain:
    def __init__(self, l, r, left_len, right_len):
        self.l = l
        self.r = r
        self.left_len = left_len
        self.right_len = right_len


def get_common_values(x_tab, y_tab):
    x_set = set()
    for col in x_tab:
        x_set.update(col)
    
    y_set = set()
    for col in y_tab:
        y_set.update(col)
    
    return x_set & y_set


def construct_hypergraph_nodes(common_values, tab):
    num_rows = len(tab[0])
    num_cols = len(tab)
    nodes = []
    col_degs = [0] * num_cols
    row_degs = [0] * num_rows
    
    # count degrees
    for c_id in range(num_cols):
        for r_id in range(num_rows):
            # only cells whose value is in common_values are kept.
            if tab[c_id][r_id] in common_values:
                col_degs[c_id] += 1
                row_degs[r_id] += 1
    
    # construct nodes                
    for c_id in range(num_cols):
        for r_id in range(num_rows):
            if tab[c_id][r_id] in common_values:
                deg = col_degs[c_id] + row_degs[r_id]
                node = HpgNode(tab[c_id][r_id], r_id, c_id, deg)
                nodes.append(node)

    return nodes


def calculate_upperbound(bidomains):
    ub = 0
    for bd in bidomains:
        ub += min(bd.left_len, bd.right_len)
    return ub


def select_bidomain(bidomains, left):
    best = -1
    best_value = float('inf')
    best_tie = float('inf')
    
    for i, bd in enumerate(bidomains):
        value = max(bd.left_len, bd.right_len)
        tie = min(left[bd.l:bd.l + bd.left_len]) if bd.left_len > 0 else float('inf')
        
        if value < best_value or (value == best_value and tie < best_tie):
            best_value = value
            best_tie = tie
            best = i

    return best


def remove_node_from_left_bidomain(left, bd, v):
    for i in range(bd.left_len):
        if left[bd.l + i] == v:
            left[bd.l + i], left[bd.l + bd.left_len - 1] = left[bd.l + bd.left_len - 1], left[bd.l + i]
            bd.left_len -= 1
            return


def remove_bidomain(bidomains, idx):
    bidomains[idx] = bidomains[-1]
    bidomains.pop()


def partition(all_indices, start, length, candidate, hg):
    i = 0
    for j in range(length):
        idx = all_indices[start + j]
        node = hg[idx]
        if node.row == candidate.row or node.col == candidate.col:
            all_indices[start + i], all_indices[start + j] = all_indices[start + j], all_indices[start + i]
            i += 1

    return i


def filter_bidomains(bidomains, left, right, x_hpg, y_hpg, v, w):
    new_bidomains = []
    for bd in bidomains:
        l = bd.l
        r = bd.r
        
        left_adj = partition(left, l, bd.left_len, v, x_hpg)
        right_adj = partition(right, r, bd.right_len, w, y_hpg)
        
        left_non_adj = bd.left_len - left_adj
        right_non_adj = bd.right_len - right_adj
        
        # only keep left and right bidomains that have the same adjacency
        if left_non_adj and right_non_adj:
            new_bidomains.append(Bidomain(l+left_adj, r+right_adj, left_non_adj, right_non_adj))
        if left_adj and right_adj:
            new_bidomains.append(Bidomain(l, r, left_adj, right_adj))
        
    return new_bidomains


def is_hyperedge_match_valid(v, w, row_mapping, col_mapping):
    if v.value != w.value:
        return False
    
    if v.row in row_mapping:
        if row_mapping[v.row] != w.row:
            return False
    else:
        if w.row in row_mapping.values():
            return False
    
    if v.col in col_mapping:
        if col_mapping[v.col] != w.col:
            return False
    else:
        if w.col in col_mapping.values():
            return False
    
    return True


def solve(x_hpg, y_hpg, best_node_maps, cur_node_maps, bidomains, left, right, row_mapping, col_mapping, delta, level, cur_level):
    # pruning
    ub = len(cur_node_maps) + calculate_upperbound(bidomains)
    best = len(best_node_maps)

    if delta is not None and level is not None and cur_level < level:
        if ub <= best * (1 + delta):
            return
    else:
        if ub <= best:
            return

    if len(cur_node_maps) > len(best_node_maps):
        best_node_maps.clear()
        best_node_maps.extend(cur_node_maps)

    bd_idx = select_bidomain(bidomains, left)
    if bd_idx == -1:
        return
    bd = bidomains[bd_idx]

    # choose a candidate node
    v_index = left[bd.l]
    v = x_hpg[v_index]
    remove_node_from_left_bidomain(left, bd, v_index)

    original_right_len = bd.right_len
    bd.right_len -= 1

    for j in range(original_right_len):
        w_index = right[bd.r + j]
        w = y_hpg[w_index]

        if not is_hyperedge_match_valid(v, w, row_mapping, col_mapping):
            continue

        new_row_mapping = row_mapping.copy()
        new_col_mapping = col_mapping.copy()

        if v.row not in new_row_mapping:
            if w.row in new_row_mapping.values():
                continue
            new_row_mapping[v.row] = w.row
        
        if v.col not in new_col_mapping:
            if w.col in new_col_mapping.values():
                continue
            new_col_mapping[v.col] = w.col

        new_bidomains = filter_bidomains(bidomains, left, right, x_hpg, y_hpg, v, w)
        cur_node_maps.append(NodePair(v_index, w_index))
        solve(x_hpg, y_hpg, best_node_maps, cur_node_maps, new_bidomains, left, right, new_row_mapping, new_col_mapping, delta, level, cur_level=cur_level+1)
        cur_node_maps.pop()

    bd.right_len += 1
    if bd.left_len == 0:
        remove_bidomain(bidomains, bd_idx)
    solve(x_hpg, y_hpg, best_node_maps, cur_node_maps, bidomains, left, right, row_mapping, col_mapping, delta, level, cur_level=cur_level+1)


@timeout(config.TIMEOUT)
def find_max_common_subhypergraph(x_hpg, y_hpg, common_values, best_node_maps, delta, level):
    left, right = [], []
    bidomains = []

    for value in common_values:
        start_left = len(left)
        start_right = len(right)
        
        for i, node in enumerate(x_hpg):
            if node.value == value:
                left.append(i)
        
        for j, node in enumerate(y_hpg):
            if node.value == value:
                right.append(j)
        
        left_len = len(left) - start_left
        right_len = len(right) - start_right
        
        if left_len > 0 and right_len > 0:
            bidomains.append(Bidomain(start_left, start_right, left_len, right_len))
    
    cur_node_maps = []
    row_mapping, col_mapping = {}, {}
    
    solve(x_hpg, y_hpg, best_node_maps, cur_node_maps, bidomains, left, right, row_mapping, col_mapping, delta, level, cur_level=0)
    return best_node_maps


def main(args):
    x_filepath = f"{args.x}.csv"
    y_filepath = f"{args.y}.csv"
    
    # load raw tables, skip specified header rows
    x_tab_df = pd.read_csv(x_filepath, skiprows=args.xh, header=None)
    y_tab_df = pd.read_csv(y_filepath, skiprows=args.yh, header=None)
    
    # convert df to lists of lists (columns), replace NaNs with empty strings
    x_tab = x_tab_df.where(pd.notna(x_tab_df), "").astype(str).values.T.tolist()
    y_tab = y_tab_df.where(pd.notna(y_tab_df), "").astype(str).values.T.tolist()

    # setup: construct hypergraphs
    setup_start_time = time.time()
    
    common_values = get_common_values(x_tab, y_tab)
    x_hpg = construct_hypergraph_nodes(common_values, x_tab)
    y_hpg = construct_hypergraph_nodes(common_values, y_tab)

    # sort nodes by degree (high to low)
    x_hpg.sort(key=lambda n: -n.deg)
    y_hpg.sort(key=lambda n: -n.deg)
    
    # solve: find maximum common subhypergraph
    solve_start_time = time.time()
    
    best_node_maps = []  # store the best node-to-node mappings found so far
    try:
        find_max_common_subhypergraph(x_hpg, y_hpg, common_values, best_node_maps, args.d, args.n)
    except TimeoutError:
        print("Timeout!")
    except Exception as e:
        print(e)
    
    end_time = time.time()

    # results
    x_row, x_col = len(x_tab[0]), len(x_tab)
    y_row, y_col = len(y_tab[0]), len(y_tab)
    overlap_size = len(best_node_maps)
    setup_time   = solve_start_time - setup_start_time
    solve_time   = end_time - solve_start_time
    total_time   = setup_time + solve_time

    print("===== Results =====")
    print(f"X Table Size   : {x_row * x_col} ({x_row} rows x {x_col} cols)")
    print(f"Y Table Size   : {y_row * y_col} ({y_row} rows x {y_col} cols)")
    print(f"Overlap Size   : {overlap_size}")
    print(f"Setup Time (s) : {setup_time:.4f}")
    print(f"Solve Time (s) : {solve_time:.4f}")
    print(f"Total Time (s) : {total_time:.4f}")


if __name__ == "__main__":
    args = parse_args()
    main(args)
