# Graph Finder
(incomplete) python wrapper for creating integer programs for finding graphs with certain constraints

## Overview

In this I use [cp_solver](https://developers.google.com/optimization/cp/cp_solver) from google's ortools in order to find strongly regular graphs with particular (n, r, λ, μ).

- `checks.py` contains short functions used to assert that the structure of a graph is as expected, like asserting that the graph is r-regular.
- `main.py` contains the solver.

## Details of particular importance

As shown below after the solver, we assert that the structure of the graph is correct, based on the emitted adjacency matrix.
```python
if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        adj = extract_adjacency_matrix(solver, vars, n)
        dump_graph(f'srg({n},{k},{lam},{mu})', adj)

        g = checks.adjacency_to_dict_rep(adj)
        assert checks.check_regular(g, k)
        assert checks.check_adjacent_regular(g, lam)
        assert checks.check_non_adjacent_regular(g, mu, n)
        render_adjacency_matrix(adj)
    else:
        print("No solution found.")
```

The only real optimization I have done is taking advantage of the adjacency matrix being symmetric.
```python
def create_graph(model, graph_size):
    vars = {}
    for i in range(graph_size):
        for j in range(i + 1, graph_size):
            vars[(i, j)] = model.NewBoolVar('x' + str(i) + str(j))
    return vars
```
