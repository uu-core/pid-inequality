from attribute_decomposition import computeInequaltyDecomposition, printDecompositition

print('-- Example 7 results ---')
f = lambda x: (x-1)**2
result = computeInequaltyDecomposition(f, 0.4, 'model.csv', ',', 'Indicator value','Number of individuals', ['Attribute 1', 'Attribute 2'])
printDecompositition(result)

print('\n-- Example 7 with GE(0.2) ---')
f = lambda x: (x**(1-0.2)-x)/(0.2*(0.2-1))
result = computeInequaltyDecomposition(f, 1e-25, 'model.csv', ',', 'Indicator value','Number of individuals', ['Attribute 1', 'Attribute 2'])
printDecompositition(result)

print('\n-- Example 7 with Atkinson index(0.8) transformation ---')
tranform = lambda x: 1-(0.8*(0.8-1)*x+1)**(1/(1-0.8))
printDecompositition([(atom, tranform(cumulative), tranform(partial)) for (atom, cumulative, partial) in result])
print('Remember the the partial contributions of a transformation are additive under a transformed addition operation (transformed inclusion-exclusion relation).')