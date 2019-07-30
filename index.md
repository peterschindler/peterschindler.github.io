---
layout: page
mathjax: false 
permalink: /
---

On this website I will try to guide you through the interactive part of the practical DFT lecture. Hopefully this will be interesting and useful for you!

We will be using the [Atomic Simulation Environment (ASE)](https://wiki.fysik.dtu.dk/ase/) in conjunction with the DFT code [GPAW](https://wiki.fysik.dtu.dk/gpaw/index.html).

### Google Colab ###
Both ASE and GPAW operate through Python. We will be using Google's Colab service to make it easy to follow along. No need to install Python 
locally - you only need a browser. Our code will be running in a Jupyter notebook on a Google server.</br>
To access a Python notebook just follow one of the links and then click on <b>'File' > 'Save a copy in Drive...'</b>
<center><img src="Images/colab.jpeg" alt="Google Colab: Copying a Python notebook to your Google Drive" style="width: 300px;"/></br>
</center>
Once you have copied the notebook to your Google Drive you can edit, execute and work on it.

###ASE Basics###
First, we need to install ASE and GPAW on the server
```bash
!apt install ase
!apt install gpaw
```
Then we can import the necessary packages into the Python notebook
```python
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from gpaw import GPAW, PW
```
This has typically one of the following forms (.module is specified only a specific module of a Python package is imported, otherwise the whole package is considered):
```python
from package.module import object/function
from package.module import *
import package.module
import package.module as pm
```
In Python we first need to import the necessary packages to


### Resources ###

1. [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase/)
2. <font color="red">color</font>