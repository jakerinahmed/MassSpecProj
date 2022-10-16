'''
MassSpecProj.py

Project idea: Allow python to read a .jdx file and plot a mass spectrum. Calculate fragments and generate a table of possible 
fragment pieces and the corresponding molecule/atom. If I have time, perhaps write my own library instead of using the JCAMP module
'''

from jcamp import JCAMP_reader #module needed to be able to read files in the JCAMP-DX format
import matplotlib.pyplot as plt #plot graphs
import numpy #perform mathematical operations on arrays and lists
from tabulate import tabulate
import periodictable as pt
import chemparse
import time

print("Welcome to MassSpecProj.py! This program is built to analyse mass spectrum data in the JCAMP-DX file format.")

while True:
    try:
        filename = input("Enter file name (Do NOT put the filename extension) : ")
        jcamp_dict = JCAMP_reader(filename + '.jdx') #A dictionary is made containing the information found in the file inputted above
        break
    except FileNotFoundError:
        print("File does not exist. Make sure you have entered the filename correctly and make sure it is saved in the same place as the python script.")
print(filename + ".jdx found")

plt.figure(figsize=(10,7)) #Plot the graph
for n in range(len(jcamp_dict['x'])):
    plt.plot((jcamp_dict['x'][n],jcamp_dict['x'][n]), 
             (0.0, jcamp_dict['y'][n]), 'm-', linewidth=1.0)
plt.title(jcamp_dict['title'])
plt.xlabel(jcamp_dict['xunits'])
plt.ylabel(jcamp_dict['yunits'])
plt.xticks(numpy.arange(min(jcamp_dict['x'])-1, max(jcamp_dict['x'])+1, 5.0))
plt.ylim(0,10000)
#x axis goes from minimum x value to maximum x value with tick spacing of 5.0

'''
A file in the .jdx format is read by the JCAMP_reader() function and creates a dictionary containing the data found in the file.
The dictionary contains the names of the headers and also two numpy arrays: x and y. These two arrays contain the x values (M/Z) and 
the y values (Relative intensity). With these two arrays, we able to plot M/Z against Relative intensity.

Source for code of plotting the graphs:
https://andthelightshattered.wordpress.com/2013/12/10/jcamp-a-python-module-for-reading-jcamp-dx-format-spectral-data/

For further information on how the jcamp module works, visit:
https://pypi.org/project/jcamp/
'''

'''
Isotopic peaks need to be removed. This has been arbitrarily set as any peak with RELATIVE INTENSITY < 200 and larger than MW.
'''

jcamp_dict['x'], jcamp_dict['y'] = zip(*((x,y) for x, y in zip(jcamp_dict['x'],jcamp_dict['y']) if y >= 200))

molecule = str(pt.formula(jcamp_dict['molform']))
MW = round(pt.formula(jcamp_dict['molform']).mass)

jcamp_dict['x'] = list(jcamp_dict['x']) #Convert tuple to list

jcamp_dict['x'] = list(filter(lambda x: x <= MW, jcamp_dict['x'])) #Only keeps values in list less than or equal to MW of molecule

plt.figure(figsize = (10,7))
for n in range(len(jcamp_dict['x'])):
    plt.plot((jcamp_dict['x'][n],jcamp_dict['x'][n]), 
             (0.0, jcamp_dict['y'][n]), 'm-', linewidth=1.0)
plt.title(jcamp_dict['title'] + ' simplified')
plt.xlabel(jcamp_dict['xunits'])
plt.ylabel(jcamp_dict['yunits'])
plt.xticks(numpy.arange(min(jcamp_dict['x'])-1, max(jcamp_dict['x'])+1, 5.0))
plt.ylim(0,10000)
#Same as code used for first graph

print("The mass of " + jcamp_dict['title'] + ' (' + molecule + ')' + " is: " + str(MW))

print("Calculating ALL possible fragments, this could take a while...")

molecule_parse = chemparse.parse_formula(molecule) #Parses the molecular formula to present the atoms and No. of atoms 

possible_fragment = None
possible_fragment_mass = 0
actual_fragments = [] #This list will contain the ACTUAL fragments
actual_fragments_masses = [] #List of masses of ACTUAL fragments

def moleculefound():
    '''
    Returns
    -------
    possiblefragment : A fragment that is randomly generated
    possiblefragmentmass : Mass of the randomly generated fragment
    '''
    random_gen_fragment_string = [] #Empty list
    
    for x,y in molecule_parse.items(): #For every key (atom) and value (number of corresponding atom) in dictionary
        y = numpy.random.randint(0,y+1) #Randomly generate the value (number of the atom) between 0 and y+1
        if y == 0:
            continue #If 0 is generated for an atom then ignore
        else:
            y = str(y) #Turn the value into a string
            random_gen_fragment_string.append(x+y) #Concatenate the atom and number of atoms and put them in a list
    
    messy_fragment = "".join(random_gen_fragment_string) #Use join function to create a "messy" chemical formula (this is a fragment)
    
    possible_fragment = pt.formula(messy_fragment) #Clean up chemical formula
    possible_fragment_mass = round(pt.formula(messy_fragment).mass) #find mass of formula

    return possible_fragment,possible_fragment_mass

t_end = 0 #Timer

if MW < 100:
    t_end = time.time() + 10
elif 100 <= MW < 200:
    t_end = time.time() + 20
elif MW >= 200:
    t_end = time.time() + 30
    
while time.time() < t_end:
    possible_fragment, possible_fragment_mass = moleculefound()
    
    for i in jcamp_dict['x']:
        molecule_already_found = False
        
        for j in actual_fragments:
            if possible_fragment == j: #If the randomly generated fragment matches a fragment that has already been found...
                molecule_already_found = True #...this boolean then becomes true
                
        if possible_fragment_mass == i and molecule_already_found == False: #If the mass of the randomly generated fragment equals to that in jcamp_dict['x'] AND an identical fragment has not yet been found...
            actual_fragments.append(possible_fragment) #...add the randomly generated fragment to the list of actual fragments
            actual_fragments_masses.append(possible_fragment_mass) #Add the mass of the actual fragment to this list to avoid having 2 fragments of the same mass
            break
            #This is to essentially, stop the same fragment appearing in the list
            
headers = ['Possible molecule','Fragment'] #Headers for the table that will be presented to user
fragment_table = zip(actual_fragments, actual_fragments_masses) #Create the table by zipping the fragments calculated before and their respective masses
fragment_table = sorted(fragment_table, key = lambda x: x[1]) #Sorts fragments by fragment mass

print(tabulate(fragment_table, headers=headers, tablefmt='fancy_grid'))