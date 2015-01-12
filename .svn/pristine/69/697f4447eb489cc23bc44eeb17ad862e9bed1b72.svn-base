#!/usr/bin/env python

import sys
from collections import OrderedDict

def join(iterable, delimeter):
    s = ""
    for i in iterable[:-1]:
        s = s + i + delimeter
    s = s + iterable[-1]
    return s

cToPy = {"int": "c_int",
          "float": "c_float",
          "double": "c_double",
          "char": "c_char",
          "const char *": "c_char_p",
          "char *": "c_char_p",
          "bool": "c_bool"
          }

cToCppPrefixes = { "string": "string(",
                    "String": "String(",
                    "char *": "strdup(",
                    "const char *": "strdup(",
                    "ResidueID": "ResidueID("
                    }
cToCppSuffixes = { "string": ")",
                    "String": ")",
                    "char *": ")",
                    "const char *": ")",
                    "ResidueID": ", ' ')"
                    }

cppToC = {"string": "const char *",
          "String": "const char *",
          "char *": "char *",
          "const char *": "const char *",
          "ResidueID": "int"
          }

wrapPrefixes = {  "string": "strdup(",
                "String": "strdup(",
                "char *": "strdup(",
                "const char *": "strdup(",
                'ResidueID':""
             }

wrapSuffixes = {   "string": ".c_str())",
                "String": ".c_str())",
                "char *": ")",
                "const char *": ")",
                "ResidueID": ".getResidueNumber()"
            }

cFileName = sys.argv[1]

structures = OrderedDict()

name = None
comment = False
for line in open(cFileName, "r"):
    line = line.strip()
    if len(line) <= 0:
        continue
    # print line
    l = line.split()
    if name and comment and (line[-1] == "}" or ")" in line):
        line = "// "+line
        structures[name][line] = ""
        comment = False
        continue
    if name and (line[-1] == "{" or "(" in line):
        comment = True
        if ")" in line:
            line = "// "+line
            structures[name][line] = ""
            comment = False
    if line.startswith("}"):
        name = None
    if line.startswith("struct") or line.startswith("class"):
        name = l[1].split("{")[0]
        if name =="MMB_EXPORT":
            name = l[2].split("{")[0]
        structures[name] = OrderedDict()
        # print "class %s(Structure):" % name
        continue
    if name and line[-1] != ";":
        line = "// "+line
        structures[name][line] = ""
        continue
    if name:
        if comment:
            line = "// "+line
            structures[name][line] = ""
            continue
        # Check if constructor/destructor
        field = line.split(";")[0]
        if len(field.split()) < 2:
            structures[name][field.split()[-1]] = "//"
            continue
        # Get type
        t = join(field.split()[:-1], " ")
        if t not in cToPy.keys()+cppToC.keys():
            print "Unknown type: " + t + " in " + name
            t = "//"+t
        # Get field name
        fieldName = field.split()[-1]
        structures[name][fieldName] = t

# C wrapper
fileOut = open(sys.argv[2], "w")      
for s in structures.keys():
    fileOut.write("typedef struct %s_wrapper{\n" % s)
    struct = structures[s]
    
    # C structure
    buff = "\tint mmbID;\n"
    for field in struct.keys():
        if field.startswith("//"):
            buff += "\t%s\n" % field
        else:
            t = struct[field]
            if t in cppToC.keys(): 
                t = cppToC[t]
            buff += "\t%s %s;\n" %(t, field)
    fileOut.write(buff + "}%s_wrapper;\n\n" % s)

    # Wrapping functions
    # void updateX_wrapper(X & _struct_, X_wrapper * _wrap_) // copy _struct_ to _wrap_
    buff = ""
    buff += "void update%s_wrapper(%s & _struct_, %s_wrapper * _wrap_){\n" % (s,s,s)
    # buff += "\t%s_wrapper _wrap_;\n" % s
    for field in struct.keys():
        if field.startswith("//"):
            buff += "\t%s\n" % field
        else:
            t = struct[field]
            if t.startswith("//"): 
                buff += "//"
            wrapPrefix = ""
            wrapSuffix = ""
            if t in cppToC: 
                wrapPrefix = wrapPrefixes.get(t,"")
                wrapSuffix = wrapSuffixes.get(t,"")
            buff += "\t_wrap_->%s = %s_struct_.%s; //%s\n" %(field, wrapPrefix, field+wrapSuffix, t)
    # buff += "\treturn _wrap_; \n}"
    fileOut.write(buff + "\n}\n")

    # void updateX(X_wrapper * _wrap_, X & _struct_) // copy _wrap_ to _struct_
    buff = ""
    buff += "void update%s(%s_wrapper * _wrap_, %s & _struct_){\n" % (s,s,s)
    for field in struct.keys():
        if field.startswith("//"):
            buff += "\t%s\n" % field
        else:
            t = struct[field]
            if t.startswith("//"): buff += "//"
            wrapPrefix = ""
            wrapSuffix = ""
            if t in cppToC: 
                wrapPrefix = cToCppPrefixes.get(t,"")
                wrapSuffix = cToCppSuffixes.get(t,"")
            buff += "\t_struct_.%s = %s_wrap_->%s; //%s\n" %(field, wrapPrefix, field+wrapSuffix, t)
    fileOut.write(buff + "\n}\n")
fileOut.close()


# Python wrapper
fileOut = open(sys.argv[3], "w")      
for s in structures.keys():
    fileOut.write("class %s_wrapper(Structure):\n" % s)
    fileOut.write("\t_fields_ = [ \n")
    struct = structures[s]
    buff = "\t\t\t('mmbID', c_int),\n"
    for field in struct.keys():
        try:
            t = struct[field]
            if t in cppToC.keys(): t = cppToC[t]
            buff += "\t\t\t('%s', %s),\n" %(field, cToPy[t])
        except KeyError:
            buff += "\t\t\t#('%s', %s),\n" %(field, struct[field]+ ": unknown ctype")
    buff = buff[:-2] + "\n"
    fileOut.write(buff + "\t\t\t]\n")
    buff = "\tdef __setattr__(self, name, value):\n"
    buff += "\t\tcmd(name + ' ' + value)\n"
    buff += "\t\tcall('updateParameterReader', self)\n"
    fileOut.write(buff+"\n")

    fileOut.write("%s_ptr = POINTER(%s_wrapper)\n" % (s, s))

    buff = "MMB.update%s_wrapper.argtypes = [c_void_p, c_char_p]\n" % s
    # buff += "MMB.update%s_wrapper.restype = c_int\n" % s
    buff += "MMB.update%s.argtypes = [c_void_p, c_char_p]\n" % s
    # buff += "MMB.update%s.restype = c_int\n" % s
    # buff += "def get%ss():\n" % s
    # buff += "\tobjs = %s_ptr()\n" % s
    # buff += "\tnbObjs = call('get%ss', byref(objs))\n" % s
    # buff += "\treturn objs[0:nbObjs]\n\n"

    fileOut.write(buff)

fileOut.close()

