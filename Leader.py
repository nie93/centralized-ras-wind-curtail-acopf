from riaps.run.comp import Component
from libopf.adapter import RscadCmdAdapter
from libopf.case import Case, Const
from libopf.compute_interface import calculatePg
from libopf.opf import runcopf
import logging
from pprint import *


__author__ = "Zhijie Nie"
__copyright__ = "Copyright (c) 2019, Zhijie Nie, Shyam Gopal, Washington State University"
__credits__ = []
__version__ = "0.1"


class Leader(Component):
    def __init__(self):
        super(Leader, self).__init__()
        self.const = Const()
        self.opf_case = Case()
        self.opf_case.import_case('./case14mod/')
        self.rscad_cmd_adapter = RscadCmdAdapter({"ip": "192.168.1.104", "port": 4575})
        self.rscad_cmd_adapter.connect_rscad()
        self.logger.info("Leader initialized")
        
    def on_clock(self):
        msg = self.clock.recv_pyobj()
        self.backuplink.send_pyobj(msg)
#         self.logger.info("Sending keepalive to backup link")
        
    def on_pmuDataReady(self):
        dframe = self.pmuDataReady.recv_pyobj()
        currentPg = calculatePg(dframe)
        self.opf_case.set_gen_prop(self.const.PMAX, [1,2,3], currentPg[1:4])
        res = runcopf(self.opf_case, flat_start=False)
                
#         self.logger.info("              Data received: %s" % str(dframe))
#         self.logger.info("      Current calculated Pg: %s" % str(currentPg))
#         self.logger.info("      Optimal results found: %s" % str(res['PG']))
        
        if res['PG'][2] < 0.99*currentPg[2]:
            self.rscad_cmd_adapter.send_cmd('SetSlider "SL3" = %.4f;' % (res['PG'][2] / 100 - 0.02))
            self.rscad_cmd_adapter.send_cmd('SetSlider "SLQG3" = %.2f;' % (res['QG'][2]))
        
    def __destroy__(self):
        self.logger.info("[%s] - Closing RSCAD connection ..." % str(self.pid))
        self.rscad_cmd_adapter.close()
        self.logger.info("[%s] - destroying RasLeader" % str(self.pid))
        