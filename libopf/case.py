import csv
import numpy as np
import os
from pdb import *


__author__ = "Zhijie Nie"
__copyright__ = "Copyright (c) 2019, Zhijie Nie, Shyam Gopal, Washington State University"
__credits__ = []
__version__ = "0.1"


class Const(object):
    def __init__(self):
        self.BUS_I       = 0    ## bus number (1 to 29997)
        self.BUS_TYPE    = 1    ## bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
        self.PD          = 2    ## Pd, real power demand (MW)
        self.QD          = 3    ## Qd, reactive power demand (MVAr)
        self.GS          = 4    ## Gs, shunt conductance (MW at V = 1.0 p.u.)
        self.BS          = 5    ## Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
        self.BUS_AREA    = 6    ## area number, 1-100
        self.VM          = 7    ## Vm, voltage magnitude (p.u.)
        self.VA          = 8    ## Va, voltage angle (degrees)
        self.BASE_KV     = 9    ## baseKV, base voltage (kV)
        self.ZONE        = 10   ## zone, loss zone (1-999)
        self.VMAX        = 11   ## maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
        self.VMIN        = 12   ## minVm, minimum voltage magnitude (p.u.)      (not in PTI format)
        self.GEN_BUS     = 0    ## bus number
        self.PG          = 1    ## Pg, real power output (MW)
        self.QG          = 2    ## Qg, reactive power output (MVAr)
        self.QMAX        = 3    ## Qmax, maximum reactive power output at Pmin (MVAr)
        self.QMIN        = 4    ## Qmin, minimum reactive power output at Pmin (MVAr)
        self.VG          = 5    ## Vg, voltage magnitude setpoint (p.u.)
        self.MBASE       = 6    ## mBase, total MVA base of this machine, defaults to baseMVA
        self.GEN_STATUS  = 7    ## status, 1 - machine in service, 0 - machine out of service
        self.PMAX        = 8    ## Pmax, maximum real power output (MW)
        self.PMIN        = 9    ## Pmin, minimum real power output (MW)
        self.PC1         = 10   ## Pc1, lower real power output of PQ capability curve (MW)
        self.PC2         = 11   ## Pc2, upper real power output of PQ capability curve (MW)
        self.QC1MIN      = 12   ## Qc1min, minimum reactive power output at Pc1 (MVAr)
        self.QC1MAX      = 13   ## Qc1max, maximum reactive power output at Pc1 (MVAr)
        self.QC2MIN      = 14   ## Qc2min, minimum reactive power output at Pc2 (MVAr)
        self.QC2MAX      = 15   ## Qc2max, maximum reactive power output at Pc2 (MVAr)
        self.RAMP_AGC    = 16   ## ramp rate for load following/AGC (MW/min)
        self.RAMP_10     = 17   ## ramp rate for 10 minute reserves (MW)
        self.RAMP_30     = 18   ## ramp rate for 30 minute reserves (MW)
        self.RAMP_Q      = 19   ## ramp rate for reactive power (2 sec timescale) (MVAr/min)
        self.APF         = 20   ## area participation factor
        self.F_BUS       = 0    ## f, from bus number
        self.T_BUS       = 1    ## t, to bus number
        self.BR_R        = 2    ## r, resistance (p.u.)
        self.BR_X        = 3    ## x, reactance (p.u.)
        self.BR_B        = 4    ## b, total line charging susceptance (p.u.)
        self.RATE_A      = 5    ## rateA, MVA rating A (long term rating)
        self.RATE_B      = 6    ## rateB, MVA rating B (short term rating)
        self.RATE_C      = 7    ## rateC, MVA rating C (emergency rating)
        self.TAP         = 8    ## ratio, transformer off nominal turns ratio
        self.SHIFT       = 9    ## angle, transformer phase shift angle (degrees)
        self.BR_STATUS   = 10   ## initial branch status, 1 - in service, 0 - out of service
        self.ANGMIN      = 11   ## minimum angle difference, angle(Vf) - angle(Vt) (degrees)
        self.ANGMAX      = 12   ## maximum angle difference, angle(Vf) - angle(Vt) (degrees)
        self.PF          = 13   ## real power injected at "from" bus end (MW)       (not in PTI format)
        self.QF          = 14   ## reactive power injected at "from" bus end (MVAr) (not in PTI format)
        self.PT          = 15   ## real power injected at "to" bus end (MW)         (not in PTI format)
        self.QT          = 16   ## reactive power injected at "to" bus end (MVAr)   (not in PTI format)        
        self.PW_LINEAR   = 0
        self.POLYNOMIAL  = 1
        self.MODEL       = 0    ## cost model, 1 = piecewise linear, 2 = polynomial 
        self.STARTUP     = 1    ## startup cost in US dollars
        self.SHUTDOWN    = 2    ## shutdown cost in US dollars
        self.NCOST       = 3    ## number breakpoints in piecewise linear cost function,
                                ## or number of coefficients in polynomial cost function
        self.COST        = 4    ## parameters defining total cost function begin in this col
                                ## (MODEL = 1) : p0, f0, p1, f1, ..., pn, fn
                                ##      where p0 < p1 < ... < pn and the cost f(p) is defined
                                ##      by the coordinates (p0,f0), (p1,f1), ..., (pn,fn) of
                                ##      the end/break-points of the piecewise linear cost
                                ## (MODEL = 2) : cn, ..., c1, c0
                                ##      n+1 coefficients of an n-th order polynomial cost fcn,
                                ##      starting with highest order, where cost is
                                ##      f(p) = cn*p^n + ... + c1*p + c0
        

class Case(object):

    def __init__(self):
        self.mva_base = 100

    def import_case(self, path):
        const = Const()

        self.path       = path
        self.bus        = np.genfromtxt(os.path.join(path, "bus.csv"), delimiter=',')
        self.gen        = np.genfromtxt(os.path.join(path, "gen.csv"), delimiter=',')
        self.branch     = np.genfromtxt(os.path.join(path, "branch.csv"), delimiter=',')
        self.gencost    = np.genfromtxt(os.path.join(path, "gencost.csv"), delimiter=',')
        self.branchrate = np.genfromtxt(os.path.join(path, "branchrate.csv"), delimiter=',')

    def set_gen_prop(self, col, idx, value):
        self.gen[idx, col] = value

    def set_branch_prop(self, type, idx, value):
        if type == 'RATE':
            self.branchrate[idx] = value
            
    def scale_branch_prop(self, col, multi):
        self.branch[:, col] = multi * self.branch[:, col]


