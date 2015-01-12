from ctypes import *
from ctypes.util import find_library

import platform

import signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

####################################################################

# Find and load the library
MMB = None
if platform.system() != "Windows":
    MMB = cdll.LoadLibrary(find_library("libMMB_Python_wrapper"))
else:
    MMB = cdll.LoadLibrary("MMB_Python_wrapper.dll")

# Retrieve MMB's functions names and pointers
functions = MMB.__dict__


####################################################################
## MB C++ exception proxy
class MMBError(Exception):
    """MMB C++ exception proxy
    
    Attributes:
        caller  -- function or parameter raising the error
        msg     -- error message from MMB
    """

    def __init__(self, caller, msg):
        self.caller = caller
        self.msg = msg

####################################################################
## Submit a command to MMB
#  See MMB's documentation for the available commands.
def cmd(cmdStr):
    """
    MMB parameter command wrapper.
    Raise a MMBError in case error
    """
    # print cmdStr
    errors = create_string_buffer("",1024)
    functions['command'](cmdStr, errors)
    if errors.value:
        msg = errors.value
        errors = create_string_buffer("",1024)
        raise MMBError(cmdStr, msg)
MMB.command.argtypes                    = [c_char_p, c_char_p]

## Call the function funcname from MMB with the arguments provided 
#  plus a string that will receive potential error messages. 
#  @param funcName the name of the function to call in string
#  @args args arguments of the function to call
#  @return the return value of the function
#  @throws MMBError in case the MMB lib raise and exception
def call(funcName, *args):
    """
    Call MMB functions by name and manage errors.
    Raise a MMBError in case of error
    """
    errors = create_string_buffer("",1024)
    args = list(args)
    args.append(errors)
    retVar = functions[funcName](*args)
    if errors.value:
        msg = errors.value
        errors = create_string_buffer("",1024)
        print funcName, msg
        raise MMBError(funcName, msg)
    return retVar

from params import *
from wrappers import *

MMBparameters = ParameterReader_wrapper()
MMBparameters.firstStage = 1
MMBparameters.lastStage = 1

