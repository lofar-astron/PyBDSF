
"""Display an animated statusbar"""
import sys
import time
from . import functions as func

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
        
        # 5 frames per second
        self.min_interval = 0.2
        self.last_print_time = 0.0
        
        self.__getsize()
        if max > 0:
            self.comp = int(float(self.pos) / self.max * self.columns)
        else:
            self.comp = 0

    # find number of columns in terminal
    def __getsize(self):
        try:
            rows, columns = func.getTerminalSize()
        except ValueError:
            rows = columns = 0
            
        if int(columns) > self.max + 2 + 44 + (len(str(self.max))*2 + 2):
            self.columns = self.max
        else:
            # note: -2 is for brackets, -44 for 'Fitting islands...' text, rest is for pos/max text
            self.columns = int(columns) - 2 - 44 - (len(str(self.max))*2 + 2)
            # Protection against negative or zero column width
            if self.columns < 1:
                self.columns = 1
        return

    # Redraw progress bar
    def __print(self):
        # Update terminal size with each frame rendering
        self.__getsize()

        sys.stdout.write('\x1b[1G')

        # Handle the case where there are no items to process (prevents from
        # dividing by zero later)
        if self.max == 0:
            sys.stdout.write(self.color + self.text + '[] 0/0\033[0m\n')
        else:
            # We are about to divide by self.max, but it is safe in the 'else' block.
            # We calculate 'self.comp' right after checking the terminal size (__getsize) 
            # so the progress bar adapts to possible window resizing.
            self.comp = int(float(self.pos) / self.max * self.columns)
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
        self.busy_char = busy_chars[self.spin_pos]

    # increment number of completed items
    def increment(self):
        self.inc = 1
        current_time = time.time()
        
        if (self.pos + self.inc) >= self.max:
            self.pos = self.max
            self.comp = self.columns
            self.busy_char = ''
            self.__print()
            return 0
        else:
            self.pos += self.inc
            self.inc = 0
            
            # Only refresh if at least 0.2 seconds have passed since the last print
            if current_time - self.last_print_time >= self.min_interval:
                self.spin()
                self.__print()
                self.last_print_time = current_time
                
        return 1

    def start(self):
        self.started = 1
        self.last_print_time = time.time()
        self.__print()

    def stop(self):
        if self.started:
            self.pos = self.max
            self.comp = self.columns
            self.busy_char = ''
            self.__print()
            sys.stdout.write('\n')
            self.started = 0
            return 0
