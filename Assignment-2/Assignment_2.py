import sys
import numpy as np
import math
import cmath

# Defining Classes for different components
class passive:
    def __init__(self, name, n1, n2, value):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = value

class volt_curr_source:
    def __init__(self, name, n1, n2, value, phase=0):
        self.name = name
        self.n1 = n1
        self.n2 = n2
        self.value = value
        self.phase = phase

# Function to parse value given in alpha numeric form to numeric form
def parse(x):
    y = len(x)
    if(not x[y-1].isalpha()):return float(x)
    if(x[y-1]=='p'):return float(x[0:y-1])* 1e-12
    if(x[y-1]=='n'):return float(x[0:y-1])* 1e-9
    if(x[y-1]=='u'):return float(x[0:y-1])* 1e-6
    if(x[y-1]=='m'):return float(x[0:y-1])* 1e-3
    if(x[y-1]=='k'):return float(x[0:y-1])* 1e3
    if(x[y-1]=='M'):return float(x[0:y-1])* 1e6
    if(x[y-1]=='G'):return float(x[0:y-1])* 1e9

if(len(sys.argv) != 2):
    print("Invalid number of arguments!")
else:
    filename = sys.argv[1]
    if(not filename.endswith(".netlist")):
        print("File type incorrect!")
    else:
        try:
            nodes = []# list for storing nodes
            node_dict = {}# dictionary for storing index of nodes
            w = 1e-100 # frequency
            SPICELines = []# list to store formatted netlist
            comp_dict = {'R': [], 'C': [], 'L': [], 'V': [],'I': []} # dictionary to store details of components
            with open(filename, "r") as f:
                #removing comments
                [SPICELines.append(line.split('#')[0].split()) for line in f.readlines()]
            if ['.circuit'] in SPICELines:
                check = 0
                for line in SPICELines:
                    if line == ['.circuit']:
                        check = 1  # If there is .circuit, start getting values
                    elif line == ['.end']:
                        check = 0  # If there is .end, stop getting values
                    elif '.ac' in line:
                        w = parse(line[2])
                    elif check:
                        # Add data to comp_dict
                        try:
                            if('R' in line[0]):
                                comp_dict['R'].append(passive(line[0], line[1], line[2], parse(line[3])))

                            elif('C' in line[0]):
                                comp_dict['C'].append(passive(line[0], line[1], line[2], parse(line[3])))

                            elif('L' in line[0]):
                                comp_dict['L'].append(passive(line[0], line[1], line[2], parse(line[3])))

                            elif('V' in line[0]):
                                if(line[3] == 'ac'):
                                    comp_dict['V'].append(volt_curr_source(line[0], line[1], line[2], parse(line[4])/(2*math.sqrt(2)), phase=parse(line[5])))
                                else:
                                    comp_dict['V'].append(volt_curr_source(line[0], line[1], line[2], parse(line[4])))

                            elif('I' in line[0]):
                                if(line[3] == 'ac'):
                                    comp_dict['I'].append(volt_curr_source(line[0], line[1], line[2], parse(line[4])/(2*math.sqrt(2)), phase=parse(line[5])))
                                else:
                                    comp_dict['I'].append(volt_curr_source(line[0], line[1], line[2], parse(line[4])))

                            if(line[1] not in nodes):
                                nodes.append(line[1])
                            if(line[2] not in nodes):
                                nodes.append(line[2])

                        except IndexError:
                            print("Value Count Incorrect!")
                            exit()

                        except ValueError:
                            print("Value Type Incorrect!")
                            exit()
                try:#append GND node to beginning of list
                    nodes.remove('GND')
                    nodes = ['GND']+nodes
                except:
                    print("No ground given!")
                    exit()

                # Creating a dictionary of index for nodes
                for i in range(len(nodes)):
                    node_dict[nodes[i]] = i

                M = np.zeros((len(nodes)+len(comp_dict['V']), len(nodes)+len(comp_dict['V'])), np.complex)
                b = np.zeros((len(nodes)+len(comp_dict['V'])), np.complex)
                # ground eqn
                M[0][0] = 1.0
                # Equations for independent voltage source
                for i in range(len(comp_dict['V'])):
                    source = comp_dict['V'][i]
                    M[len(nodes)+i][node_dict[source.n1]] = -1.0
                    M[len(nodes)+i][node_dict[source.n2]] = 1.0
                    b[len(nodes)+i] = cmath.rect(source.value,source.phase*cmath.pi/180)
                # Voltage source equations
                if(len(comp_dict['V']) != 0):
                    for v in comp_dict['V']:
                        if(node_dict[v.n1] != 0):
                            M[node_dict[v.n1]][len(nodes)] = -1.0
                        if(node_dict[v.n2] != 0):
                            M[node_dict[v.n2]][len(nodes)] = 1.0
                # Resistor equations
                if(len(comp_dict['R']) != 0):
                    for r in comp_dict['R']:
                        if(node_dict[r.n1] != 0):
                            M[node_dict[r.n1]][node_dict[r.n1]] += 1/r.value
                            M[node_dict[r.n1]][node_dict[r.n2]] -= 1/r.value
                        if(node_dict[r.n2] != 0):
                            M[node_dict[r.n2]][node_dict[r.n1]] -= 1/r.value
                            M[node_dict[r.n2]][node_dict[r.n2]] += 1/r.value
                # Capacitor equations
                if(len(comp_dict['C']) != 0):
                    for c in comp_dict['C']:
                        if(node_dict[c.n1] != 0):
                            M[node_dict[c.n1]][node_dict[c.n1]] += complex(0, 2*np.pi*w*c.value)
                            M[node_dict[c.n1]][node_dict[c.n2]] -= complex(0, 2*np.pi*w*c.value)
                        if(node_dict[c.n2] != 0):
                            M[node_dict[c.n2]][node_dict[c.n2]] += complex(0, 2*np.pi*w*c.value)
                            M[node_dict[c.n2]][node_dict[c.n1]] -= complex(0, 2*np.pi*w*c.value)
                # Inductor equations
                if(len(comp_dict['L']) != 0):
                    for l in comp_dict['L']:
                        if(node_dict[l.n1] != 0):
                            M[node_dict[l.n1]][node_dict[l.n1]] -= complex(0, 1/(w*2*np.pi*l.value))
                            M[node_dict[l.n1]][node_dict[l.n2]] += complex(0, 1/(w*2*np.pi*l.value))
                        if(node_dict[l.n2] != 0):
                            M[node_dict[l.n2]][node_dict[l.n2]] -= complex(0, 1/(w*2*np.pi*l.value))
                            M[node_dict[l.n2]][node_dict[l.n1]] += complex(0, 1/(w*2*np.pi*l.value))
                # equations for Independent Current source
                if(len(comp_dict['I']) != 0):
                    for source in comp_dict['I']:
                        if(node_dict[source.n1] != 0):
                            b[node_dict[source.n1]] += cmath.rect(source.value, source.phase*cmath.pi/180)
                        if(node_dict[source.n2] != 0):
                            b[node_dict[source.n2]] -= cmath.rect(source.value, source.phase*cmath.pi/180)
                x = np.linalg.solve(M, b)
                i=0
                for n in nodes:
                    print(" Node Voltage " + n+" =  " , x[i])
                    i = i+1
                for v in comp_dict['V']:
                    print("Current Through Voltage Source " + v.name + " = ", x[i])
                    i= i+1
            else:
                print('Netlist not compatible to given format!')
        except FileNotFoundError:
            print("Given file does not exist")
