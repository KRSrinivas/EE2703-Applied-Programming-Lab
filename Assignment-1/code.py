import sys

# Extracting the parameter from branch
def branch_parameter(branch):
    parameter = branch.split()

    # R, L, C, Independent Sources
    if(len(parameter) == 4):
        element = parameter[0]
        node1 = parameter[1]
        node2 = parameter[2]
        value = parameter[3]
        return [element, node1, node2, value]

    # Current Controlled -X- Source
    elif(len(parameter) == 5):
        element = parameter[0]
        node1 = parameter[1]
        node2 = parameter[2]
        voltageSource = parameter[3]
        value = parameter[4]
        return [element, node1, node2, voltageSource, value]

    # Voltage Controlled -X- Source
    elif(len(parameter) == 6):
        element = parameter[0]
        node1 = parameter[1]
        node2 = parameter[2]
        voltageSourceNode1 = parameter[3]
        voltageSourceNode2 = parameter[4]
        value = parameter[5]
        return [element, node1, node2, voltageSourceNode1, voltageSourceNode2, value]

    else:
        return []

#function to print netlist as given
def printcktstrt(Branchparameters):
    for x in Branchparameters[:]:
        for y in x[:]:
            print(y, end=' ')
        print('')
    print('')
    return

#function to print the netlist in reverse
def printcktrev(Branchparameters):
    for x in Branchparameters[::-1]:
        for y in x[::-1]:
            print(y, end=' ')
        print('')
    print('')
    return

# checking number of command line arguments
if len(sys.argv)!=2 :
    sys.exit("Invalid number of arguments!")
else:
    try:
        circuitFile = sys.argv[1]

        # checking if given netlist file is of correct type
        if (not circuitFile.endswith(".netlist")):
            print("File type incorrect!")
        else:
            with open (circuitFile, "r") as f:
                SPICELines = []
                for line in f.readlines():
                    SPICELines.append(line.split('#')[0].split('\n')[0])
                try:
                    # finding the location of the identifiers
                    start = SPICELines.index(".circuit")
                    end = SPICELines.index(".end")

                    Branches = SPICELines[start+1:end]
                    Branchparameters = [branch_parameter(line) for line in Branches]

                    # Printing Circuit Definition in Reverse Order
                    print("\nThe Circuit Netlist is:\n")
                    printcktstrt(Branchparameters)

                    print("\nThe Reversed Circuit Netlist is:\n")
                    printcktrev(Branchparameters)

                except ValueError:
                    print("Netlist not compatible to given format!")
    except FileNotFoundError:
        print("Given file does not exist")
