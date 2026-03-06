# SPICE — Spatial Parametric Inference for Correspondence Estimation

**SPICE (Spatial Parametric Inference for Correspondence Estimation)** presents an analytical approach that infers the statistical significance of spatial association tests between spatial maps by estimating the **effective degrees of freedom** of the data.

---

## Installation

We provide both MATLAB and Python implementations of SPICE.

### MATLAB

Download the `matlab_spice` folder from the present repository and add it to your MATLAB path:

```matlab
addpath(genpath('matlab_spice'))
```

### Python

The Python implementation is distributed via PyPI.  
You may also clone and install manually from the present repository.

> **Note**  
> The package is distributed on PyPI under the name **`brainspice`**,  
> but imported in Python as **`spice`** after installation.


#### PyPI installation

```bash
pip install brainspice
```

#### Import in Python

```python
import spice
```

---

## Example usage

To evaluate the association between spatial maps *x* and *y* given:

- spatial map ***x***, may contain NaN or Inf values.
- spatial map ***y***, may contain NaN or Inf values.
- spatial coordinates ***coord*** of observations

### SPICE (stationary assumption)

#### MATLAB

```matlab
% take 70–80s to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
% pef - significance p-value
% rX - Pearson correlation coefficient
% nef - effective sample size
[pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2] = effective_sample_size_estimation(x, y, coord);
```

#### Python

```python
# take ~3.5 min to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
# pef - significance p-value
# rX - Pearson correlation coefficient
# nef - effective sample size

import spice
pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2 = spice.effective_sample_size_estimation(x, y, coord)
```

### SPICE-NS (nonstationary assumption)

#### MATLAB

```matlab
% SPICE-NS with data-driven parcellation
% take 70–80s to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map
[pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2] = effective_sample_size_estimation(x, y, coord, 'xparc', 'auto', 'yparc', 'auto');
```

#### Python

```python
# SPICE-NS with data-driven parcellation
# take ~3.5 min to run (Apple Silicon M1 Pro) on fsaverage5 10k cortical map

import spice
pef, rX, nef, run_status, n_parc, p_naive, fc_para1, fc_para2 = spice.effective_sample_size_estimation(x, y, coord,xparc='auto',yparc='auto')
```

---

## Documentation

For full documentation and additional examples, see
[SPICE documentation](https://brainspice.readthedocs.io/en/latest/index.html)
