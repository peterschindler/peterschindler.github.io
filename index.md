---
layout: page
mathjax: false 
permalink: /
---

On this website I will try to guide you through the interactive part of the practical DFT lecture. Hopefully this will be interesting and useful for you!

We will be using the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/) in conjunction with the DFT code [GPAW](https://wiki.fysik.dtu.dk/gpaw/index.html).

## Google Colab ##
Both ASE and GPAW operate through Python. We will be using Google's Colab service to make it easy to follow along. No need to install Python 
locally - you only need a browser. Our code will be running in a Jupyter notebook on a Google server.<br/>
To access a Python notebook just follow [**this link**](https://colab.research.google.com/drive/1mofb8yD9rcwVFlBnwkrMf2rx3E1Pru6u) and then click on <b>'File' > 'Save a copy in Drive...'</b>
<center><img src="Images/colab.jpeg" alt="Google Colab: Copying a Python notebook to your Google Drive" style="width: 300px;"/><br/>
</center>
Once you have copied the notebook to your Google Drive you can edit, execute and work on it.

## ASE Basics ##
### Importing modules ###
First, we need to install ASE and GPAW on the server (`%%capture` represses output of the notebook cell)
```bash
%%capture
!apt install ase
!apt install gpaw
```
Then we can import the necessary packages into the Python notebook
```python
from ase import Atoms #The Atoms object is used to define and work with atomic structure in ASE 
from ase.build import bulk #The ase.build module has many functions to generate solids and molecules - very handy!
from ase.io import read, write #The ase.io module is used to read and write crystal/molecular structure files
from gpaw import GPAW, PW #GPAW will be our main DFT calculator and PW is the plane wave mode
```
This has typically one of the following forms (if .module is specified only a specific module of a Python package is imported, otherwise the whole package is considered):
```python
import package.module as pm #let's you access functions/objects like so: pm.function
from package.module import object/function #imports object or function from package.module
from package.module import * #imports all functions/objects
import package.module #same
```
This is sorted from best to worst practice when it comes to importing modules in Python.
### Creating our first atomic structure ###
First, let's create a H2 molecule (experimental distance is 0.74 Angstroms) and add some vacuum padding around it
```python
h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.74)])
h2.center(vacuum=2.5)
```
To export the atomic structure we created we can use the ase.io module
```python
write('h2.cif', h2)
```
and if we wanted to we can let ASE tell us the cell parameters and atomic positions of the object we created
```python
h2.cell
h2.positions
```
### Set GPAW calculator ###
Now, let's define a calculator using GPAW. We are specifying the exchange-correlation functional, *k*-points, the plane wave energy cutoff and a file to safe the calculation results.
```python
calc = GPAW(xc='LDA', 
            kpts=(1,1,1), 
            mode=PW(500), 
            txt='h2.txt')
```
Then we need to "attach" our calculator to the H2 molecule we created earlier.
```python
h2.set_calculator(calc)
```
### Evaluate properties ###
Now we can simply calculate the energy of the atomic object and could directly print it or safe the result in a variable
```python
h2.get_potential_energy()
```
alternatively, we could also get the forces by `.get_forces()` among many other things.
### Testing for convergence ###
Python's for loop can now we very useful in testing for convergence (without tediously rewriting input files or continuously reading out energy results)
```python
cutoffs = [100,200,300,400,500,1000]
energies = []

for cutoff in cutoffs:
  calc = GPAW(xc='LDA', 
              kpts=(1,1,1), 
              mode=PW(cutoff), 
              txt='h2.txt')
  h2.set_calculator(calc)
  energies.append(h2.get_potential_energy())
print(energies)#eV
```
here `cutoffs` is a list of energy cutoffs and the for loop iterates over each entry and then adds it to the (at first empty) `energies` list.
### Resources ###

1. [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
2. <font color="red">color</font>