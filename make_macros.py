import sys

config_file = open("Config")

names = []
values = []
types = []

#
# Get the configuration options
# and figure out what is a float
#

for line in config_file:
    if (line[0]=='#'):
        pass
    else:
        tokens = line.split('=')
        if (len(tokens) < 2):
            pass
        else:
            isFloat = 0
            
            names.append(tokens[0])
            #if len(tokens > 1):
            values.append(tokens[1])
            #else:
            #    values.append("NULL")
            
            for char in tokens[1]:
                if (char == '.'):
                    isFloat = 1
            if (isFloat == 1):
                types.append(1)
            else:
                types.append(0)

#print values
#print names
#print types

config_file.close()
header = open("./include/hllc_defs.h",'w')
i_max = len(values)
i=0

while (i < i_max):
    if (types[i] == 0 and values[i] != "NULL"):
        header.write("#define " + str(names[i]) + " " + str(int(values[i])) + '\n')
    elif(types[i] != 0 and values[i] != "NULL"):
        header.write("#define " + str(names[i]) + " " + str(float(values[i])) + '\n')
    else:
        header.write("#define " + str(names[i]))
    i = i + 1
