<img style="float:right;width:100px;" src="../img/scarecrow.png" alt="scarecrow"/>

[Back to root](../README.md)

# Installation

To install `scarecrow` first create a conda/mamba environment. After creating the environment, activate it, and install the dependencies.

```bash
mamba create --name scarecrow python=3.12
mamba activate scarecrow
mamba install git numpy pandas pip psutil pympler pysam rich scipy seaborn
pip install git+https://github.com/MorganResearchLab/scarecrow.git
pip install pyahocorasick
```

Whenever using `scarecrow` ensure that the environment is first activated (`mamba activate scarecrow`).
