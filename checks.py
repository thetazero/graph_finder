def adjacency_to_dict_rep(A):
    g = {}
    for i in range(len(A)):
        for j in range(len(A)):
            if A[i][j] == 1:
                if i not in g:
                    g[i] = set()
                if j not in g:
                    g[j] = set()
                g[i].add(j)
                g[j].add(i)
    return g

def check_regular(g, r):
    """
    Checks if the graph g is r-regular. 
    """
    for i in g:
        if len(g[i]) != r:
            print(g[i])
            return False
    return True

def check_adjacent_regular(g, lam):
    """
    Check if every pair of adjacent vertices have lam common neighbors. 
    """
    for i in g:
        for j in g[i]:
            common = set(g[i]).intersection(set(g[j]))
            if len(common) != lam:
                return False
    return True

def check_non_adjacent_regular(g, lam, n):
    """
    Check if every pair of adjacent vertices have lam common neighbors. 
    """
    for i in g:
        all = set(range(n))
        for j in g[i]:
            all.remove(j)
        all.remove(i)
        for j in all:
            common = set(g[i]).intersection(set(g[j]))
            if len(common) != lam:
                return False
    return True
