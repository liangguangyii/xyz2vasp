## VASP to xyz and xyz to VASP

### xyz to vasp
```python
file_name="Au9Pt4+"          #* do NOT include .xyz
eleNames = ['Au', 'Pt']
boxadd = 20
#! in INCAR the order of elements is very important!!!

xyz2vasp(file_name, eleNames, boxadd)
```

### vasp to xyz
```python
vasp2xyz('CONTCAR')     #* file name must be CONTCAR or POSCAR
```