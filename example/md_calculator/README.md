# Example: Minima Hopping with two calculators
In some cases a numerically cheaper second calculator for the molecular dynamics part and a pre-optimization might be useful. In this example presented here two calculators are initialized. One is attached tot the atoms object and the other is given as an argument to the MinimaHopping class. The MD and pre-optimization are performed by using latter while the data base is built with the calculator attached to the atoms object. The geometry optimization in this example is performed using the Quantum Espresso calculator. The MD and a pre-optimization is performed using openKIM.

More information on how to use the above described calculators can be found in the ASE documentation (https://wiki.fysik.dtu.dk/ase/ase/calculators/calculators.html). The openKIM perfoms the MD using a Stillinger-Weber potential (SW_StillingerWeber_1985_Si__MO_405512056662_005) and can be installed via
```bash
kim-api-collections-management install user SW_StillingerWeber_1985_Si__MO_405512056662_005
```
after having installed the openKIM API correctly. Quantum Espresso is described in detail on https://www.quantum-espresso.org/. It is possible to installed via sudo:
```bash
sudo apt install quantum-espresso
```
However, if you want to calculate more complex systems than the example described here performance might be crucial. 



The example can be started by executing:
```bash
python mh_na35_second_calculator.py
```
