# pid-inequality

This repository provides an implementation of the method in:
- Tobias Mages and Christian Rohner. _Quantifying Redundancies and Synergies with Measures of Inequality_

## Overview
- `demo.py`: provides a simple usage example
- `attribute_decomposition.py`: implements the presented method

_**Requirements:**_ [pandas](https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html)

## Interface
`attribute_decomposition.py`: `computeInequaltyDecomposition(.)` take the followign arguments:
- `f`: the generator function for an f-divergence
- `p`: the parameter `0<=p<=1`
- `modelCSV`: the CSV-file containing the system model
- `delimiter`: the used deliminator in the CSV-file
- `indicatorColumnName`: the column name  in the CSV-file for the indicator variable.
- `individualCountColumnName`: the column name in the CSV-file for the number of individuals with this indicator and attributes.
- `AttributeColumnNameList`: the list of column name in the CSV-file for the attributes of the decomposition.
- `return value`: a list of tuples. Each tuple contains the atom name, cumultive and partial result.

## Example usage
```
from attribute_decomposition import computeInequaltyDecomposition, printDecompositition

# define the generator function of an f-divergence
f = lambda x: (x-1)**2

# compute the decomposition
result = computeInequaltyDecomposition(f, 0.4, 'model.csv', ',', 'Indicator value','Number of individuals', ['Attribute 1', 'Attribute 2'])

# print the result
printDecompositition(result)
#               Atom               , cumulative, partial
# [['Attribute 1', 'Attribute 2']] ,   0.589246, 0.000000
#        [['Attribute 2']]         ,   0.299376, 0.117162
#        [['Attribute 1']]         ,   0.360000, 0.056538
#[['Attribute 1'], ['Attribute 2']],   0.416538, 0.172708
#               [[]]               ,   0.000000, 0.242838
```
