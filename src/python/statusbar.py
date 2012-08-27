"""Display an animated statusbar"""
import sys
import os
import functions as func

class StatusBar():
    # class variables:
    # max:  number of total items to be completed
    # pos:  number of completed items
    # spin_pos: current position in array of busy_chars
    # inc:  amount of items to increment completed 'pos' by
    #           (shared resource)
    # comp: amount of '=' to display in the progress bar
    # started: whether or not the statusbar has been started
    # color: color of text 
    def __init__(self, text, pos=0, max=100, color='\033[0m'):
        self.text = text
        self.pos = pos
        self.max = max
        self.busy_char = '|'
        self.spin_pos = 0
        self.inc =  0
        self.started = 0
        self.color = color
        self.__getsize()
        if max > 0:
            self.comp = int(float(self.pos) / self.max * self.columns)
        else:
            self.comp = 0
            
    # find number of columns in terminal
    def __getsize(self):
        try:
#             rows, columns = os.popen('stty size', 'r').read().split()
            rows, columns = func.getTerminalSize()
        except ValueError:
            rows = columns = 0
        if int(columns) > self.max + 2 + 44 + (len(str(self.max))*2 + 2):
            self.columns = self.max
        else:
            # note: -2 is for brackets, -44 for 'Fitting islands...' text, rest is for pos/max text
            self.columns = int(columns) - 2 - 44 - (len(str(self.max))*2 + 2)
        return

    # redraw progress bar
    def __print(self):
        self.__getsize()

        sys.stdout.write('\x1b[1G')
        if self.max == 0:
            sys.stdout.write(self.color + self.text + '[] 0/0\033[0m\n')
        else:
            sys.stdout.write(self.color + self.text + '[' + '=' * self.comp + self.busy_char + '-'*(self.columns - self.comp - 1) + '] ' + str(self.pos) + '/' + str(self.max) + '\033[0m')
            sys.stdout.write('\x1b[' + str(self.comp + 2 + 44) + 'G')
        sys.stdout.flush()
        return

    # spin the spinner by one increment
    def spin(self):
        busy_chars = ['|','/','-','\\']
        self.spin_pos += 1
        if self.spin_pos >= len(busy_chars):
            self.spin_pos = 0
        # display the busy spinning icon
        self.busy_char = busy_chars[self.spin_pos]
        sys.stdout.write(self.color + busy_chars[self.spin_pos] + '\x1b[1D' + '\033[0m')
        sys.stdout.flush()
        
    # increment number of completed items
    def increment(self):
        self.inc = 1
        if (self.pos + self.inc) >= self.max:
            self.pos = self.max
            self.comp = self.columns
            self.busy_char = ''
            self.__print()
            sys.stdout.write('\n')
            return 0
        else:
            self.pos += self.inc
            self.inc = 0
            self.spin()
            self.comp = int(float(self.pos) / self.max \
                * self.columns)
            self.__print()
        return 1

    def start(self):
        self.started = 1
        self.__print()

