from ortools.sat.python import cp_model
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from timeit import default_timer as timer


def order_small_big(x):
    (a, b) = x
    if a < b:
        return (a, b)
    return (b, a)


def r_regularity_constraint(model, vars, n, r):
    for i in range(n):
        edges = []
        for j in range(n):
            if i != j:
                edges.append(vars[order_small_big((i, j))])
        model.Add(sum(edges) == r)


def get_edges(vars, from_vertex, n):
    edges = []
    for j in range(n):
        if from_vertex != j:
            edges.append(vars[order_small_big((from_vertex, j))])
    return edges

def and_constraint(model, a, b, c):
    """
    Adds a constraint that a AND b == c
    
    Trick From:
    https://github.com/google/or-tools/blob/stable/ortools/sat/docs/boolean_logic.md#product-of-two-boolean-variables
    """
    model.AddBoolOr(a.Not(), b.Not(), c)
    model.AddImplication(c, a)
    model.AddImplication(c, b)

def make_common_neighbor_vars(model, vars, n):
    """
    Creates a new variable l(i, j, k) which is equal to 
    l <=> i ~ k AND j ~ k
    In other words, l is true iff i and j have a common neighbor k.
    """
    new_vars = {}
    for (i, j) in vars:
        for k in range(n):
            if k not in [i, j]:
                l = model.NewBoolVar(f'l{i}{j}{k}') 
                new_vars[(i, j, k)] = l

                ik_edge = vars[order_small_big((i, k))]
                jk_edge = vars[order_small_big((j, k))]

                and_constraint(model, ik_edge, jk_edge, l) # l <=> i_edges[k] AND j_edges[k]
    return new_vars


def adjacent_regularity_constraint(model, vars, common_neighbor_vars, n, lam):
    """
    Every two adjacent vertices have lambda common neighbors.
    """
    for (i, j) in vars:
        # l = i_edges[k] * j_edges[k] trick via
        # l <=> i_edges[k] AND j_edges[k]
        ij_edge = vars[order_small_big((i, j))]

        shared_edges = []
        for k in range(n):
            if k not in [i, j]:
                l = common_neighbor_vars[(i, j, k)]

                shared_edges.append(l)

        common = sum(shared_edges)

        model.Add(common == lam).OnlyEnforceIf(ij_edge)

def extract_adjacency_matrix(solver, vars, n):
    adj = np.zeros((n, n))
    for (i, j) in vars:
        adj[i][j] = solver.Value(vars[(i, j)])
        adj[j][i] = solver.Value(vars[(i, j)])

    return np.array(adj)


def render_adjacency_matrix(adjacency_matrix):
    gr = nx.from_numpy_array(adjacency_matrix)
    nx.draw(gr, node_size=500, with_labels=True)
    plt.show()


def create_graph(model, graph_size):
    vars = {}
    for i in range(graph_size):
        for j in range(i + 1, graph_size):
            vars[(i, j)] = model.NewBoolVar('x' + str(i) + str(j))
    return vars


def example():
    model = cp_model.CpModel()
    solver = cp_model.CpSolver()

    n = 10  # Number of vertices
    k = 3  # Every vertex has k neighbors
    lam = 0 # Every two adjacent vertices have lam common neighbors
    mu = 1  # Every two non-adjacent vertices have mu common neighbors
    vars = create_graph(model, n)

    r_regularity_constraint(model, vars, n, k)
    common_neighbor_vars = make_common_neighbor_vars(model, vars, n)
    adjacent_regularity_constraint(model, vars, common_neighbor_vars, n, lam)

    solver.parameters.cp_model_probing_level=0

    start = timer()
    status = solver.Solve(model)
    end = timer()
    print(end-start)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adj = extract_adjacency_matrix(solver, vars, n)
        # for (i, j, k) in l_vars:
        #     print(f"l({i},{j},{k}) = {solver.Value(l_vars[(i, j, k)])}")
        render_adjacency_matrix(adj)
    else:
        print("No solution found.")

if __name__ == '__main__':
    example()
