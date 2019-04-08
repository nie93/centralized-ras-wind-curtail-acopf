from riaps.run.comp import Component
import logging
import os
import threading
import time
from pypmu import pdc
import struct
from apt_offline_core.AptOfflineMagicLib import NONE


__author__ = "Zhijie Nie"
__copyright__ = "Copyright (c) 2019, Zhijie Nie, Shyam Gopal, Washington State University"
__credits__ = []
__version__ = "0.1"


class C37118InputDataAdapter(Component):
    
    def __init__(self):
        super().__init__()
        self.pid = os.getpid()
        self.conn_dict = [{"ip": "192.168.1.111", "port": 4712, "id": 1},
                          {"ip": "192.168.1.111", "port": 4722, "id": 2}]
#         self.conn_dict = [{"ip": "192.168.1.132", "port": 4712, "id": 12},
#                           {"ip": "192.168.1.132", "port": 4722, "id": 12},
#                           {"ip": "192.168.1.119", "port": 4712, "id":  9}]
#         self.conn_dict = [{"ip": "192.168.1.119", "port": 4722, "id":  9}]
        self.logger.info("[%s] - initialized C37118InputDataAdapter" % str(self.pid))
        self.c37118Thread = None
        
    def on_clock(self):
        now = self.clock.recv_pyobj()
        self.logger.info("Sending keepalive: %s" % now)
        if self.c37118Thread == None:
            self.c37118Thread = C37118Thread(self)
            self.c37118Thread.start()
            self.queue.activate()
        
    def __destroy__(self):
        self.logger.info("[%s] - destroying C37118InputDataAdapter" % str(self.pid))
        
    def on_queue(self):
        dataFrame = self.queue.recv_pyobj()
        self.logger.info("[%s] - queueing C37118InputDataAdapter" % str(self.pid))
        self.pmuData.send_pyobj(dataFrame)


class C37118Thread(threading.Thread):
    def __init__(self, component):
        threading.Thread.__init__(self)
        self.port = component.queue
        self.active = threading.Event()
        self.active.clear()
        self.waiting = threading.Event()
        self.terminated = threading.Event()
        self.terminated.clear()
        self.component = component
        self.pdc = list()
        self.phnmr = dict()
        self.chnme = dict()
        self.dframes = dict()
        self.data = dict()

    def run(self):
        self.plug = self.port.setupPlug(self)             # Ask parent port to make a plug for this end

        if self.terminated.is_set(): return

        self.active.wait()

        for conn in self.component.conn_dict:
            pmu = pdc.Pdc(pmu_ip=conn['ip'], pmu_port=conn['port'], pdc_id=conn['id'])
            try:
                pmu.run()
                self.pdc.append(pmu)
            except:
                pass
            
        for pmu in self.pdc:
            self.parse_config(pmu)
            pmu.start()

        while 1:
            if self.terminated.is_set():
                for pmu in self.pdc:
                    pmu.quit()
                break
 
            if not self.active.is_set():
                for pmu in self.pdc:
                    pmu.stop()
                    self.active.wait()
                    pmu.start()
 
            self.get_dframes()
            phasors = self.parse_dframes()
            self.plug.send_pyobj(phasors)     
            
    def get_dframes(self):
        if not self.pdc:
            return None
        
        _dframes = dict()
        for pmu in self.pdc:
            pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
            _dframes[pmukey] = pmu.get()
        
        self.dframes = _dframes
        
        return

    def parse_config(self, pmu):
        if not pmu:
            return None

        pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
        self.chnme[pmukey] = list()
        config = pmu.get_config()
        phnmr = struct.unpack('>h', config[40:42])[0]
        self.phnmr[pmukey] = phnmr

        chnam = ""
        for iph in range(phnmr):
            chnam = ""
            for i in range(46 + 16 * iph, 46 + 16 * (iph+1)):
                chnam += struct.unpack('>c', config[i:i+1])[0].decode("utf-8")
            self.chnme[pmukey].append(chnam.replace(" ", ""))

        return

    def parse_dframes(self):
        if not self.dframes:
            return None

        _dframes = self.dframes
        data = dict()
        for pmu in self.pdc:
            pmukey = "%s:%d" % (pmu.pmu_ip, pmu.pmu_port)
            data[pmukey] = dict()
            
            for iph in range(self.phnmr[pmukey]):
                magkey = self.chnme[pmukey][iph] + "_mag"
                angkey = self.chnme[pmukey][iph] + "_ang"
                data[pmukey][magkey] = struct.unpack('>f', _dframes[pmukey][16 + 8*iph:20 + 8*iph])[0]
                data[pmukey][angkey] = struct.unpack('>f', _dframes[pmukey][20 + 8*iph:24 + 8*iph])[0]
                      
        return data

    def activate(self):
        self.active.set()

    def deactivate(self):
        self.active.clear()

    def terminate(self):
        self.terminated.set()

    