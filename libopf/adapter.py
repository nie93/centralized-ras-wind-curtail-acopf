import struct
import socket


__author__ = "Zhijie Nie"
__copyright__ = "Copyright (c) 2019, Zhijie Nie, Shyam Gopal, Washington State University"
__credits__ = []
__version__ = "0.1"


class RscadCmdAdapter(object):

    def __init__(self, connection_dict):
        self.s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.conn_dict = connection_dict
        
    def connect_rscad(self):
        self.s.connect((self.conn_dict['ip'], self.conn_dict['port'])) # RSCAD socket location
        self.send_cmd('Start;')
        self.send_cmd('ListenOnPortHandshake("temp_string");')

        # Echo to synchronize the adapter with RSCAD
        rmsg = self.s.recv(64).decode("utf-8")
        while ('temp_string' not in rmsg):
            rmsg = self.s.recv(64).decode("utf-8")
        print("*    Successfully connected to RSCAD Script TCP Server.")

    def send_cmd(self, cmdstr):
        print("*    Sending command: %s" % cmdstr)
        self.s.send(cmdstr.encode())

    def close(self):
        self.s.shutdown(socket.SHUT_RDWR)
        self.s.close()
        print("*    Successfully closed RSCAD Script TCP Server connection.")


