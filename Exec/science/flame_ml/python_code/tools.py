import sys
import torch

class Logger(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.terminal = sys.stdout
        self.log = open(log_file, "a")
        self.log.write("MaestroFlame\n")
        self.log.close()

    def write(self, message):
        self.terminal.write(message + '\n')
        self.log = open(self.log_file, "a")
        self.log.write(message)
        self.log.close()

    def flush(self):
        pass
