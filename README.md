# JuliaMolSim_Docs

These docs show how to fit and regularise PoSH potentials in Julia, and also show compatibility through ASE in Python. After installing Julia one should add the JuliaMolSim repository.

Form the Julia REPL:

`] registry add https://github.com/JuliaMolSim/MolSim.git`

From there we can add the dependencies:

`] add PoSH`

`] add IPFitting`

`] add JuLIP`

For Python support:

First we will need to install PyJulia

`pip install pyjulia`

Which requires PyCall (in Julia) to be built with your Python execubtable. Afterwards get `git clone` and `python setup.py install` this repo:

`https://github.com/casv2/pyjulip`




