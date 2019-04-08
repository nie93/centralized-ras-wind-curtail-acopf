# [RIAPS-APPS-WSU] AcopfWindCurtailmentApp - Washington State University

This application is a **Remedial Action Scheme (RAS) application**, that minimizes the wind generation curtailment when the power grid is subject to operational violations.

## Introduction

This application formulates the minimal wind generation curtailment problem in the modified IEEE 14-Bus Test System (extended with three wind farms), with AC powerflow equality constraints, together with voltage boundaries, and MVA line ratings inequality constraints. This application acquires synchrophasor data packets under IEEE C37.118 protocols from Real-Time Digital Simulator (RTDS) GTNET-PMU card. The simulation is interfaced with RSCAD running a listener script based on socket connection, that accomodates the RAS control action determined by this application.

## Developers

* Zhijie Nie [nie@ieee.org](mailto:nie@ieee.org)
* Shyam Gopal [shyam.gopal@wsu.edu](mailto:shyam.gopal@wsu.edu)

## Package Requirements

* Python 3.x
* `numpy==1.16.1`, `scipy==0.19.1`
* `pyPMU` ((Included in the repository) [https://github.com/iicsys/pypmu](https://github.com/iicsys/pypmu)

## Hardware Requirements

* Real-Time Digital Simulator (RTDS)
* GTNET-PMU card


