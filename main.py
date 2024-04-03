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

def adjacent_regularity_constraint(model, vars, n, lam):
    """
    Every two adjacent vertices have lambda common neighbors.
    """

def extract_adjacency_matrix(solver, vars, n):
    adj = [[0 for i in range(n)] for j in range(n)]
    for (i,j) in vars:
        adj[i][j] = solver.Value(vars[(i,j)])
        adj[j][i] = solver.Value(vars[(i,j)])
    
    return np.array(adj)

def render_adjacency_matrix(adjacency_matrix):
    gr = nx.from_numpy_array(adjacency_matrix)
    nx.draw(gr, node_size=500)
    plt.show()

def create_graph(model, graph_size):
    vars = {}
    for i in range(graph_size):
        for j in range(i + 1, graph_size):
            vars[(i, j)] = model.NewIntVar(0, 1, 'x' + str(i) + str(j))
    return vars

def example():
    model = cp_model.CpModel()
    solver = cp_model.CpSolver()

    n = 99 # Number of vertices
    r = 14 # Every vertex has 14 neighbors
    vars = create_graph(model, n)
    
    r_regularity_constraint(model, vars, n, r)

    start = timer()
    status = solver.Solve(model)
    end = timer()
    print(end-start)

    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adj = extract_adjacency_matrix(solver, vars, n)
        render_adjacency_matrix(adj)
    else:
        print("No solution found.")

if __name__ == '__main__':
    example()
