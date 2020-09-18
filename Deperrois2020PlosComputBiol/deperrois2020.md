
Short-term depression and long-term plasticity together tune sensitive range of synaptic plasticity
==============================

Find here code and information related to

**Deperrois N and Graupner M (2019).**
Short-term depression and long-term plasticity together tune sensitive range of synaptic plasticity.
*accepted in PLoS Comput Biol*; [bioRxiv 565291; doi: 10.1101/565291](https://doi.org/10.1101/565291).
* [Find related scripts and description here](Deperrois2019PlosComputBiol/deperrois2019.md)



Figures
-----------

### J Neurosci 2016 Fig. 1, 3 and 5 : Plasticity for irregular spike pairs.
<img src="outputFigures/Graupner2016JNeurosci_linearCaModel.png" width="600px" />

### J Neurosci 2016 Fig. 9 : Plasticity for irregular spike pairs in the calcium based model with nonlinear calcium dynamics.
<img src="outputFigures/Graupner2016JNeurosci_nonlinearCaModel.png" width="600px" />


Change parameters and re-produce Figures
-----------

Equivalently, for the calcium dynamics model with short-term plasticity for the somatosensory cortex plasticity data (Markram et al. 1997, Science)
```python
python irregularPairs_stp-Ca-Makram-Model.py
python Graupner2018_STPMarkramCaModel.py
```

Or for the calcium dynamics model with short-term plasticity for the visual cortex plasticity data (Sjoestroem 2001, Neuron)
```python
python irregularPairs_stp-Ca-Sjoestroem-Model.py
python Graupner2018_STPSjoestroemCaModel.py
```


