# Dissipative particle dynamics to reproduce conformation

This implementation of DPD is based on the repository (https://github.com/KPavelI/dpd) which was used for the simulation chromatin interactions with the nuclear lamina [1].

## Hardaware requirements

The RAM, HDD and CPU are depend on the system size and available simulation time. Current calculations were performed on cluster lomonosov-2 [2], for each run we used 56 cores.

## Software requirements

Fortran90 compiler. The software has been tested on the following systems: `Ubuntu 16.04`, `CentOS Linux` (release 7.1.1503).

## General notes

To visualize output restart files you may convert it to **mol2** using `rst2mol2.py` from the mentioned [repository](https://github.com/KPavelI/dpd) from example folder. If you are using this code, please, cite the folowing works [3, 4].

## References.

1. Ulianov, S. V., Doronin, S. A., Khrameeva, E. E., Kos, P. I., Luzhin, A. V., Starikov, S. S., ... & Mikhaleva, E. A. (2019). Nuclear lamina integrity is required for proper spatial organization of chromatin in Drosophila. Nature communications, 10(1), 1176.
2. Voevodin, V. V., Antonov, A. S., Nikitenko, D. A., Shvets, P. A., Sobolev, S. I., Sidorov, I. Y., ... & Zhumatiy, S. A. (2019). Supercomputer Lomonosov-2: large scale, deep monitoring and fine analytics for the user community. Supercomputing Frontiers and Innovations, 6(2), 4-11.
3. Gavrilov, A. A.; Chertovich, A. V.; Khalatur, P. G.; Khokhlov, A. R. Effect of Nanotube Size on the Mechanical Properties of Elastomeric Composites. Soft Matter 2013, 9 (15), 4067–4072.
4. Groot, R. D., & Warren, P. B. (1997). Dissipative particle dynamics: Bridging the gap between atomistic and mesoscopic simulation. The Journal of chemical physics, 107(11), 4423-4435.
