# Example: Minima Hopping with two calculators
In some cases a numerically cheaper second calculator for the molecular dynamics part and a pre-optimization might be useful. In this example presented here two calculators are initialized. One is attached tot the atoms object and the other is given as an argument to the MinimaHopping class. The MD and pre-optimization are performed by using latter while the data base is built with the calculator attached to the atoms object. In the example presented here the Lennard-Jones calcualtor is the second calculator while the EAM calculator is attached to the atoms object. 
The example can be started by executing:
```bash
python mh_na35_second_calculator.py
```
Note that this is only an example to show how this code feature works. More realistic examples are a DFT calculator with looser settings or a machine learned potential to perform the MD and pre-optimization. Another example are also biased runs where the second calculator includes a bias whereas the geometry optimization occurs on the PES without bias.

