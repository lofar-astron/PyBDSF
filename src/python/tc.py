"""Defines some basic facilities for handling typed values.


It's quite basic and limited implementation tailored specifically for
use in the PyBDSM user-options and derived properties. For a user
option, one can define a group that is used when listing the options to
the screen. For a property (e.g., flux density), one can define the
column name to be used on output and the associated units.

For a much more generic and capable implementation I can recommend
to look at Enthought Traits package: 
    http://code.enthought.com/projects/traits


Defined are:
 - a number tc-handlers which allow to type-check and/or cast
   values to the specific type (tcCType, tcEnum, tcTuple, 
   tcOption, tcInstance, tcList, tcAny). These aren't really
   inteded for use by end-user.

 - class TC, which implements a concept of type-checked property
   with default value.

 - a number of wrappers around TC to simplify it's usage (Int,
   Float, Bool, String, Tuple, Enum, Option, NArray, Instance,
   tInstance, List, Any)

Usage:
For the most needs it's enough to use wrapper-interface. 
One important remark -- class containing tc-variables should be
new-style class, thus you should explicitly inherit from 'object'
for Python < 2.6.

Example:
from tc import Int, Float, Bool, String, Tuple, Enum, \\
    Option, NArray, Instance, Any, TCInit

class tst(object):
    intval = Int(doc="Integer value")
    boolval = Bool(True, "Some boolean flag")
    op_type = Enum("op1", "op2", doc="Some enumerated value")

    def __init__(self):
        TCInit(self)    ### this is optional

v = tst()
v.intval = 1           # OK
v.intval = "33"        # OK, casted to 33
v.intval = "failure"   # FAILS
v.op_type= "op2"       # OK
v.op_type= "op3"       # FAILS
"""
import exceptions
import types

_sequence_types = (types.ListType, types.TupleType)
_class_types = (types.ClassType, types.TypeType)
_basic_types = (types.BooleanType, types.IntType, types.LongType,
                types.FloatType,   types.ComplexType,
                types.StringType,  types.UnicodeType)


############################################################
## Wrappers around TC to simplify it's usage for end-users
############################################################
def Int(value=0, doc=None, group=None, colname=None, units=None):
    """Create tc-value of type int"""
    return TC(value, tcCType(int), doc, group, colname, units)

def Float(value=0., doc=None, group=None, colname=None, units=None):
    """Create tc-value of type float"""
    return TC(value, tcCType(float), doc, group, colname, units)

def Bool(value=False, doc=None, group=None):
    """Create tc-value of type bool"""
    return TC(value, tcCType(bool), doc, group)

def String(value='', doc=None, group=None, colname=None, units=None):
    """Create tc-value of type string"""
    return TC(value, tcCType(str), doc, group, colname, units)

def Tuple(*values, **kws):
    """Create tc-value of type tuple.

    Parameters:
    values: zero or more arguments
    kws: keyword arguments. Currently only 'doc' and 'group'
    are recognized

    If the first item of values is a tuple, it's used as the
    default value. The remaining arguments are used to build
    type constraints and should be TC values.

    Examples:
    Tuple((1,2,3))          # tuple of 3 integers, default = (1,2,3)
    Tuple(Int(3), Float(2)) # tuple of int&float, default = (3, 2.0)
    Tuple((1,2), Int(3), Float(2)) # tuple of int+float, default = (1, 2.0)
    """
    doc = kws.pop('doc', None)
    group = kws.pop('group', None)
    if len(values) == 0:
        return TC((), tcTuple(), doc, group)

    default = None
    if isinstance(values[0], tuple):
        default, values = values[0], values[1:]

    if default is None:
        default = tuple([x._default for x in values])

    if len(values) == 0:
        values = [tc_from(x) for x in default]

    return TC(default, tcTuple(*values), doc, group)

def Enum(*values, **kws):
    """Create tc-value of type enum.

    Parameters:
    values: list or tuple of valid values
    kws: keyword arguments. Currently only 'doc'  and 'group'
    are recognized

    Default value is taken to be values[0].

    Examples:
    Enum(3, [1,2,3])    # enum of 1,2,3 with default of 3
    Enum(1,2,3)         # enum of 1,2,3 with default of 1
    """
    default = values[0]
    if (len(values) == 2) and (type(values[1]) in _sequence_types):
        values = values[1]

    doc = kws.pop('doc', None)
    group = kws.pop('group', None)

    return TC(default, tcEnum(*values), doc, group)

def Option(value, type=None, doc=None, group=None):
    """Creates optional tc-value.

    Parameters:
    value, type: default value and type
    doc: doc-string for the value
    group: group designation for the value
    """
    if type is None:
        type = tc_from(value)

    if isinstance(value, TC):
        value = value._default

    return TC(value, tcOption(type), doc, group)

def NArray(value=None, or_none=True, doc=None, group=None, colname=None, 
           units=None):
    """Creates tc-value which holds Numpy arrays

    Parameters:
    value: default value
    or_none: if 'None' is valid value
    group: group designation for the value
    colname: name of column if quantity is to be output
    units: units if quantity is to be output
    """
    try:
        import numpy as N
    except:
        raise tcError, "Can't create tc-value of type NArray " \
            "without access to numpy module"

    return Instance(value, N.ndarray, or_none, doc, group, colname, units)

def Instance(value, type=None, or_none=True, doc=None, group=None, 
             colname=None, units=None):
    """Creates tc-value which holds instances of specific class.

    Parameters:
    value, type: default value and type
    or_none: flag if 'None' is valid value for this variable
    group: group designation for the value
    colname: name of column if quantity is to be output
    units: units if quantity is to be output

    Examples:
    Instance(instance, class)
    Instance(instance)
    Instance(class)
    """
    if type is None:
        if isinstance(value, _class_types):
            value, type = None, value
        else:
            type = value.__class__

    return TC(value, tcInstance(type, or_none), doc, group, colname, units)

def tInstance(type, or_none=False):
    """Create tc-handler for values which are instances of
    the specific class.

    This function is useless on it's own, and should be
    used to create Instane-constrain for compound tc-values.
    It's especially usefull for classes which have non-trivial
    constructors.

    Parameters:
    type: target type/class
    or_none: flag if 'None' is valid value for this variable

    Example: we want to define tc-variable holding a list of objects
    List(Instance(slice, or_none=False) ## FAILS, no default value
    List(Instance(slice))  ## works, but list MAY contain None's
    List(tInstance(slice)) ## GOOD
    """
    if not isinstance(type, _class_types):
        type = type.__class__

    return tcInstance(type, or_none)

def List(value, type=None, doc=None, group=None, colname=None, units=None):
    """Creates tc-value which represents a list, where each element
    obeys specific type-constrains.

    Parameters:
    doc: docstring for the object
    value, type: default value and type
    group: parameter group to which the option belongs
    colname: name of column if quantity is to be output
    units: units if quantity is to be output


    Examples:
    List(Int())        # list of integers, default value is []
    List([1,2], Int()) # list of integers, default value is [1,2]


    Just one more warning -- List always has default value
    ([] in the simples case), and this default value is shared
    between the instances, so be carefull to not modify it.

    Counter-example for it:
    class tst(object):
        l = List(Int())

    x1 = tst()
    x2 = tst()   # both instances share default value

    x1.l.append(1)
    print x2.l   # this will print [1]

    x1.l = [2]
    print x2.l   # still [1], as x1 has it's own local value now    
    """
    if type is None:
        value, type = [], tc_from(value)

    return TC(value, tcList(type), doc, group, colname, units)

def Any(value=None, doc=None, group=None):
    """Creates tc-value of arbitrary type
    (e.g. no type-checking is done)
    """
    return TC(value, tcAny(), doc, group)

def TCInit(obj):
    """Initialize tc-variables in the new instance"""
    TC.set_property_names(obj.__class__)
    obj._tc_values = {}


############################################################
## Exception type
############################################################
class tcError(exceptions.Exception):
    """Custom exception type to simplify exception handling"""
    pass


############################################################
## TC -- type-checked variable
############################################################
class TC(object):
    """TC is an implementation of the typed-checked value.

    The primary usage pattern is via class attributes:

    class Test(object):   ### MUST be new-style object
        value1 = Int(3)
        value2 = Tuple(Int(5), Option(Any()))

    test = Test()
    print test.value1
    test.value2 = (3, None)

    An important restriction -- it might only be used with
    new-style objects (e.g. objects derived from 'object' 
    or 'type'. And the attribute should be defined in the
    class of the object.
    """
    def __init__(self, value, _type=None, doc=None, group=None, colname=None,
                 units=None):
        """Create typed-checked object.

        Parameters:
        value: default value
        _type: type specification (instance of tcHandler) or None
        doc: docstring for the object
        group: parameter group to which the option belongs
        colname: name of column if quantity is to be output
        units: units if quantity is to be output
        """
        if _type is not None:
            self._type = _type
        else:
            self._type = tc_from(value)

        self._default = self._type.cast(value)
        self._name = None # name is unknown atm
        self._group = group
        self._doc = doc
        self._colname = colname
        self._units = units

        self.__doc__ = "default value is %s (%s)" % \
                (str(self._default), self._type.info())

        if doc is not None:
            self.__doc__ += "\n" + doc

    def __get__(self, instance, cls):
        """Get a value from instance (or return default value)"""
        if instance is None:
            return self

        try:
            return instance._tc_values[self]
        except:
            return self._default

    def __set__(self, instance, value):
        """Set a value"""
        try:
            values = instance._tc_values
        except:
            values = instance._tc_values = {}

        if not self._name:
            self.set_property_names(instance.__class__)
        
        values[self] = self._type.cast(value, self._name, 
                                       instance.__class__.__name__)

    def __delete__(self, instance):
        """Revert value to default"""
        try:
            del instance._tc_values[self]
        except:
            pass

    def cast(self, value, *args):
        """See tcHandler.cast"""
        return self._type.cast(value, *args)

    def info(self):
        """Return description of tc-value"""
        return self.__doc__

    def doc(self):
        """Return short description of tc-value"""
        return self._doc

    def group(self):
        """Return group designation of tc-value"""
        return self._group

    def colname(self):
        """Return column name designation of tc-value"""
        return self._colname

    def units(self):
        """Return units designation of tc-value"""
        return self._units

    @staticmethod
    def set_property_names(klass):
        """Scan class definition and update _name for all
        TC objects defined there"""
        for k,v in klass.__dict__.iteritems():
            if isinstance(v, TC):
                v._name = k


############################################################
## tcHandler and derived handlers for the specific 
## types/values
############################################################
class tcHandler(object):
    """Base class for all tc-handlers"""
    def cast(self, value, *args):
        """Check that provided value meets type requirements
        or cast it to the specific type.
        """
        self.error(strx(value), *args)

    def is_valid(self, value):
        """Check if provided value can be safely casted to the 
        proper type"""
        try:
            self.cast(value)
            return True
        except:
            return False

    def info(self):
        """A description of a valid values"""
        return "value of unknown type"

    def error(self, value, *args):
        if len(args) == 2 and args[0]:
            error = "Failed to set property %s of class %s " \
                "to a value of %s; expected %s." % \
                (args[0], args[1], value, self.info())
        else:
            error = "A value of %s can't be casted to %s" % \
                (value, self.info())
        raise tcError(error, value, self.info(), *args)


############################################################
class tcAny(tcHandler):
    """Allows any values of any type"""
    def cast(self, value, *args):
        return value

    def info(self):
        return "any value"


############################################################
class tcCType(tcHandler):
    """Ensures that value has a specific python type

    This handler implements so-called casting-approach, where
    it will accept all values which can be converted to the
    required type by the means of casting operation. For
    example:

    v = tcCType(int)
    print v.cast(3)      # casted to 3
    print v.cast(3.3)    # casted to 3
    print v.cast("3")    # casted to 3
    """
    def __init__(self, _type):
        """Creates tcType handler.
        
        Parameters:
        _type: Python type object or a value of a reqired type
        """
        if not isinstance(_type, types.TypeType):
            _type = type(_type)

        self.type = _type

    def cast(self, value, *args):
        if type(value) is self.type:
            return value

        try:
            return self.type(value)
        except:
            self.error("%s (%s)" % (str_type(value), reprx(value)),
                       *args)

    def info(self):
        return "a value of %s" % str_type(self.type)
                    

############################################################
class tcEnum(tcHandler):
    """Ensures that a value is a member of a specified list of values"""
    def __init__(self, *values):
        """Creates a tcEnum handler.
        
        Parameters:
        values: list or tuple of all legal values

        Description:
        The list of values can be provided as a list/tuple of values
        or just specified in-line. So that ''tcEnum([1,2,3])'' and
        ''tcEnum(1,2,3)'' are equivalent.
        """
        if len(values) == 1 and type(values[0]) in _sequence_types:
            values = values[0]

        self.values = values

    def cast(self, value, *args):
        if value in self.values:
            return value

        self.error(repr(value), *args)

    def info(self):
        res = "a value of %s" % \
            " or ".join([repr(x) for x in self.values])
        return res


############################################################
class tcTuple(tcHandler):
    """Ensures that a value is a tuple of specified length,
    with elements that are of specified type
    """
    def __init__(self, *args):
        """Creates a tcTuple handler.

        Parameters:
        args: list of tuple components

        Description:
        Each tuple component should be either a specific 
        tc-handler or a value which can be converted to it
        (by the means of tc_from function)
        """
        self.tcs = tuple([tc_from(x) for x in args])

    def cast(self, value, *args):
        try:
            if type(value) in _sequence_types:
                if len(value) == len(self.tcs):
                    res = []
                    for i, h in enumerate(self.tcs):
                        res.append(h.cast(value[i]))
                    return tuple(res)
        except:
            pass

        self.error(reprx(value), *args)

    def info(self):
        res = "a tuple of the form: (%s)" % \
            ", ".join([x.info() for x in self.tcs])
        return res


############################################################
class tcOption(tcHandler):
    """Implements an optional value: None or a value 
    restricted by another tcHandler"""
    def __init__(self, _type):
        """Creates tcOption handler.
        
        Parameters:
        _type: tc-handle, Python type object or a value of 
               a reqired type
        """
        self.type = tc_from(_type)

    def cast(self, value, *args):
        try:
            if value is None:
                return value
            return self.type.cast(value)
        except:
            self.error("%s (%s)" % (str_type(value), reprx(value)),
                       *args)

    def info(self):
        return self.type.info() + " or None"


############################################################
class tcInstance(tcHandler):
    """Ensures that a value belongs to a specified python 
    class or type (or one of it's subclasses).
    """
    def __init__(self, klass, or_none=True):
        """Creates tcInstance handler.

        Parameters:
        klass: Python class, type or an instance of python class
        or_none: whether we should accept None as a valid value
                 (defaults to True)
        """
        if not isinstance(klass, _class_types):
            klass = klass.__class__
        self.klass = klass
        self.or_none = or_none

    def cast(self, value, *args):
        if (value is None) and self.or_none:
            return value
        if isinstance(value, self.klass):
            return value

        self.error(reprx(value), *args)

    def info(self):
        res = "an instance of " + str_type(self.klass)
        if self.or_none:
            res += " or None"

        return res


############################################################
class tcList(tcHandler):
    """Ensures that a value is a list containing elements of
    a specified kind. It also ensures that any change made
    to the list does't violate the list type constrains.
    """
    def __init__(self, kind):
        """Creates tcList handler.

        Parameters:
        kind: tc-handler constraining elements of the list
        """
        self.type = tc_from(kind)

    def cast(self, value, *args):
        if isinstance(value, _sequence_types):
            return tcListObject(self, value, args)

        self.error(reprx(value), *args)

    def info(self):
        return "a list where each element is " + self.type.info()


############################################################
class tcListObject(list):
    """Helper class for tcList.

    It's basically a customized implementation of list type,
    which imposes specific type constrains on it's elements.
    """
    def __init__(self, tc, values, extras):
        self.list_handler = tc
        self.type = tc.type
        self.extras = extras

        ## type-check initial values
        self.__setslice__(0, 0, values)

    def __setitem__(self, key, value):
        v = self.type.cast(value, *self.extras)
        list.__setitem__(self, key, v)

    def __setslice__(self, i, j, values):
        cast = self.type.cast
        v = [cast(x, *self.extras) for x in values]
        list.__setslice__(self, i, j, v)

    def append(self, value):
        v = self.type.cast(value, *self.extras)
        list.append(self, v)

    def extend(self, values):
        cast = self.type.cast
        v = [cast(x, *self.extras) for x in values]
        list.extend(self, v)

    def insert(self, idx, value):
        v = self.type.cast(value, *self.extras)
        list.insert(self, idx, v)


############################################################
def tc_from(v):
    """tc_from tries to guess an appropriate tc-handler for the
    provided object. 

    The basic logic is a following:
     - TC object results in it's internal type constrain
     - for a instances and type-objects of the basic numerica 
       types we use tcCType handler
     - a list of values results in tcEnum handler
     - a tuple of values results in tcTuple handler
     - a value of None results in tcAny handler
    """
    if isinstance(v, TC):
        return v._type
    if isinstance(v, tcHandler):
        return v
    if v in _basic_types:
        return tcCType(v)
    if type(v) in _basic_types:
        return tcCType(v)
    if type(v) is types.ListType:
        return tcEnum(v)
    if type(v) is types.TupleType:
        return tcTuple(*v)
    if v is None:
        return tcAny()

    error = "Can't create tc-handler for a value of %s (%s)" %\
        (str_type(v), reprx(v))
    raise tcError(error)


############################################################
def str_type(v):
    """Pretty-print type of v"""
    if isinstance(v, _class_types):
        return repr(v)[1:-1]
    else:
        return repr(type(v))[1:-1]


############################################################
def reprx(v):
    """Pretty-print value of v"""
    if type(v) is types.InstanceType:
        return v.__class__.__name__
    else:
        return repr(v)
