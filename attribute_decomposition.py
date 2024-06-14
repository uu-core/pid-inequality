import pandas as pd
from functools import reduce, cache
from itertools import product, chain, combinations

# Note: sotchastic matrices are represented as list of tuples (x,y)

# helper function: extract stochastic matrix from model
def stochasticMatrixFromModel(modelCSV,delimiter,indicatorColumnName,individualCountColumnName, AttributeColumnNameList):
    AttributeColumnNameList = list(set(AttributeColumnNameList)) # remove duplicates
    model = pd.read_csv(modelCSV,delimiter=delimiter, engine='python')
    if AttributeColumnNameList == []:
        return [(1,1)]
    else:
        partition = model.groupby(AttributeColumnNameList)
        def aggregate_subgroups(subgroup):
            groupsize = subgroup[individualCountColumnName].sum()
            indicator = (subgroup[indicatorColumnName] * subgroup[individualCountColumnName]).sum()
            return pd.Series({indicatorColumnName: indicator/groupsize, individualCountColumnName: groupsize})
        population = partition.apply(aggregate_subgroups, include_groups=False).reset_index()
        size = sum(population[individualCountColumnName])
        indSum = sum([row[individualCountColumnName]*row[indicatorColumnName] for (_, row) in population.iterrows()])
        return [(row[individualCountColumnName]/size, row[individualCountColumnName]*row[indicatorColumnName]/indSum) for (_, row) in population.iterrows()]

# helper function: compute zonogon-join of two matrices
def convexHull(matrix1, matrix2):
    # clip functions in case of numerical inaccuracy
    clipMax = lambda x: (min([1, x[0]]), min([1, x[1]])) 
    clipMin = lambda x: (max([0, x[0]]), max([0, x[1]])) 
    def integ(matrix): # sum sorted vector segments to zonogon points
        acc = (0,0)
        return [acc := clipMax((acc[0]+p[0],acc[1]+p[1])) for p in matrix]
    def diff(curve):   # diff zonogon points to vector segments
        return [curve[0]] + [clipMin((curve[i][0]-curve[i-1][0],curve[i][1]-curve[i-1][1])) for i in range(1,len(curve))]
    zono1 = integ(sorted(matrix1,key=lambda x: x[1]/x[0] if x[0] != 0 else float("inf")))
    zono2 = integ(sorted(matrix2,key=lambda x: x[1]/x[0] if x[0] != 0 else float("inf")))
    slope = lambda v1, v2: (v2[1]-v1[1])/(v2[0]-v1[0])
    def comvexHullx(acc, current, other):
        while current != []: # walk along x, always picking next point with minimal slope since both points are sorted, we can pick greedy
            if acc[-1][0] >= other[0][0]: # the next point of the other curve is behind (drop it)
                other = other[1:]
            elif slope(acc[-1], current[0]) < slope(acc[-1], other[0]): # pick the next point with minimal slope to current position
                acc, current, other = acc + [current[0]], current[1:], other
            else:
                acc, current, other = acc + [other[0]], other[1:], current
        return acc
    return diff(comvexHullx([zono1[0]], zono1[1:], zono2) if slope((0,0), zono1[0]) < slope((0,0), zono2[0]) else comvexHullx([zono2[0]], zono2[1:], zono1))

# helper function: compute join for a list of matrices
def bigJoin(matrixList):
    return reduce(convexHull, matrixList, [(1,1)])

# helper function: f-inequality
def fineq(f,p,matrix):
    r = lambda v: (p*v[0] + (1 - p)*v[1])*f(v[0]/(p*v[0] + (1 - p)*v[1]))
    return sum([r(v) for v in matrix if v != (0,0)])

# helper function: f-inequality cumulative measure
def fineqCup(stochasticMatrixGen, f, p, atom):
    return fineq(f, p, bigJoin([stochasticMatrixGen(attributeList) for attributeList in atom]))

'''
Description: compute the decomposition for a given f-inequality and model.
- `f`: the generator function for an f-divergence
- `p`: the parameter `0<=p<=1`
- `modelCSV`: the CSV-file containing the system model
- `delimiter`: the used deliminator in the CSV-file
- `indicatorColumnName`: the column name  in the CSV-file for the indicator variable.
- `individualCountColumnName`: the column name in the CSV-file for the number of individuals with this indicator and attributes.
- `AttributeColumnNameList`: the list of column name in the CSV-file for the attributes of the decomposition.
- `return value`: a list of tuples. Each tuple contains the atom name, cumultive and partial result.
'''
def computeInequaltyDecomposition(f,p,modelCSV,delimiter,indicatorColumnName,individualCountColumnName, AttributeColumnNameList):
    assert 0 <= p and p <= 1
    if p == 0:
        print('Warning: substituded p=0 with p=1e-20. The use of p=0 may require simplifying the function r(f,p) based on the specific function f to avoid a division by zero.')
        p = 1e-20
    stochasticMatrixGen = lambda x: stochasticMatrixFromModel(modelCSV, delimiter, indicatorColumnName, individualCountColumnName, x)
    reduceAtom = lambda cmp, atom: [list(set(a)) for a in atom if (not (any([cmp(a,b) for b in atom])))]
    complement = lambda atom: [[b for b in AttributeColumnNameList if b not in a] for a in atom]
    dual = lambda atom: [list(x) for x in reduceAtom(lambda x,y: set(x) > set(y), [list(set(a)) for a in product(*complement(atom))])]
    powerset = lambda attr: [list(x) for x in chain.from_iterable(combinations(attr, r) for r in range(len(attr)+1))]
    def powerSetFiltered(xs):
        if xs == []:
            return [[]]
        else:
            ys = powerSetFiltered(xs[1:])
            return [t for t in ys + [[xs[0]] + y for y in ys] if all([not (set(s1) < set(s2)) for s1 in t for s2 in t])]
    atoms = [x for x in powerSetFiltered(powerset(AttributeColumnNameList)) if x != []] # construct all lattice atoms

    @cache
    def fineqCupCached(atom): # memoize the cumulative measure (taking frozensets)
        return fineqCup(stochasticMatrixGen, f, p, [list(x) for x in atom])
    fineqCupCachedList = lambda atom : fineqCupCached(frozenset([frozenset(x) for x in atom])) # wrapper for type conversion

    res = []
    for atom in atoms:
        cumulative = fineqCupCachedList(atom)
        partial = sum([(-1)**(len(beta)-1)*fineqCupCachedList(reduceAtom(lambda x,y: set(x) < set(y), beta+atom)) for beta in powerset(dual(atom))]) if atom != [AttributeColumnNameList] else 0
        partial = round(partial,15) + 0 # round to decent precision and add zero to remove pythons 'negative zero'. 
        res.append((atom, cumulative, partial))
    return res

'''
Description: print the decomposition results
- `resultTuples`: a list of tuples. Each tuple contains the atom name, cumultive and partial result.
- `return value`: None
'''
def printDecompositition(resultTuples):
    maxLen = max([len(str(x[0])) for x in resultTuples])
    print(f"{'Atom'.center(maxLen,' ')}, {'cumulative'}, {'partial'}")
    for (atom, cumulative, partial) in resultTuples:
        print(f"{str(atom).center(maxLen,' ')},   {cumulative:.6f}, {partial:.6f}")
