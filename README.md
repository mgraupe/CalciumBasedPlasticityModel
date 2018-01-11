
Calcium-based plasticity model
==============================

In the calcium-based plasticity model, pre- and postsynaptic spike induce calcium transients. Synaptic depression or potentiation is induced whenever the compound calcium trace crosses the depression or the potentiation threshold, respectively. 

The here published python scripts implement the calculations to obtain the change in synaptic strength for pre-post spike-pairs, pre-spike and post-pair, as well as 
irregular pair stimulation. 

For more details, please refer to :

**Graupner M and Brunel N (2012).**
Calcium-based plasticity model explains sensitivity of synaptic changes to spike pattern, rate, and dendritic location. 
[*PNAS 109 (10): 3991-3996.*](http://www.pnas.org/content/109/10/3991.abstract)

**Graupner M, Wallisch P and Ostojic S (2016).**
Natural Firing Patterns Imply Low Sensitivity of Synaptic Plasticity to Spike Timing Compared with Firing Rate. 
[*J Neurosci 36(44):11238-11258*](http://www.jneurosci.org/content/36/44/11238)


Features
-----------
* The `timeAboveThreshold` class calculates the time the compound calcium trace spends above a given threshold for pre-post spike-pair, pre-spike and post-burst, as well as irregular spike pairs stimulation protocols. It furthermore calculates the time above threshold for the non-linear calcium model with regular and irregular spike pair stimulation (see Graupner et al. 2016, Fig. 9).  
* The up and down transition probabilities are furthermore calculated and converted into a change in synaptic strength considering the initial distribution of synapses and the ratio of synaptic strength between the UP and the DOWN state. 
* The basic results are plotted. 

### PNAS 2012 Fig. 2 : Diversity of STDP curves in response to spike pair stimulation.
<img src="outputFigures/Graupner2012PNAS_Fig2.png" width="400px" />

### PNAS 2012  Fig. 3 : Numbers of postsynaptic spikes and repetitions of the stimulation motif qualitatively change the STDP curve.
<img src="outputFigures/Graupner2012PNAS_Fig3.png" width="200px" />

### PNAS 2012 Fig. 4 : Plasticity for spike pairs vs. firing frequency.
<img src="outputFigures/Graupner2012PNAS_Fig4B.png" width="300px" />

### J Neurosci 2016 Fig. 1, 3 and 5 : Plasticity for irregular spike pairs.
<img src="outputFigures/Graupner2016JNeurosci_linearCaModel.png" width="600px" />

### J Neurosci 2016 Fig. 9 : Plasticity for irregular spike pairs in the calcium based model with nonlinear calcium dynamics.
<img src="outputFigures/Graupner2016JNeurosci_nonlinearCaModel.png" width="600px" />

Change parameters and re-produce Figures
-----------

**PNAS 2012:**
The parameter values for the stimulation protocol can be changed in `Graupner2012PNAS_FigXX.py` (with XX=2,3 or 4B). Equivalently, the parameters of the plasticity model implementation can be changed in  `synapticChange.py` . The png and pdf versions of the figures can be produced by running the scripts
```python
python Graupner2012PNAS_Fig2.py
python Graupner2012PNAS_Fig3.py
python Graupner2012PNAS_Fig4B.py
```

**J Neurosci 2016:** The data has to be first calculated by running `irregularPairs_Ca-Model.py` for the linear calcium dynamics model and `irregularPairs_nonlinear-Ca-Model.py` for the nonlinear calcium dynamics model. For the linear calcium dynamics model, the calcium integral is calculated semi-analytically in a C++ program, which has to be compiled first by running `make` in the `timeAboveThreshold` directory. The compound figures are subsequently generated by the scripts `Graupner2016JNeurosci_linearCaModel.py` and `Graupner2016JNeurosci_nonlinearCaModel.py`. In other words, the png and pdf versions of the linear calcium dynamics figures can be produced by running the following two scripts
```python
python irregularPairs_Ca-Model.py
python Graupner2016JNeurosci_linearCaModel.py
```

Equivalently, for the nonlinear calcium dynamics model
```python
python irregularPairs_nonlinear-Ca-Model.py
python Graupner2016JNeurosci_nonlinearCaModel.py
```

Equivalently, for the calcium dynamics model with short-term plasticity
```python
python irregularPairs_stp-Ca-Model.py
python irregularPairs_STPCaModel.py
```

Requires
-----------
Standard python packages such as **numpy**, **scipy**, **pylab**, **time**, **os**,  **sys** and **matplotlib** are required.

License
-----------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

