# -*- coding: latin-1 -*-
"""
Name: Risa Hocking
File: RamanReader.py
------------------
This class takes in some path to a Raman or PL file (map or single spot) and
returns the data at each point in a large dictionary. Will auto-detect and sort
the units into energy or wavenumbers as well as calculate the other set of values.

Coding is in latin-1 due to the use of some greek characters in the files
from the Horiba Raman tools (Xplora, Labram).
"""


import numpy as np
import os


class RamanMap():
	
    def __init__(self, file):
        self.f = file # filepath
        self.fname = self.get_filename()
        self.wn = None # ndarray of wavenumbers, shape (N, 1)
        self.energy = None # ndarray of energy, shape (N, 1)
        self.position = None # ndarray of tuples of (X, Y) in µm for each scan, shape (M,1)
        self.intensity = None # ndarray of peak intensities for each positions, shape (M x N)
        self.import_data()


    def get_filename(self):
        """
        Takes in the entire path to a file. Returns just
        the filename as a string (no preceeding part of the path).
        """
        return (os.path.split(self.f)[1])[:-4]


    def tsvtoarr(self, string, dataset):
        """
        Converts a string with tab-separated values into 
        a list. Appends this new list into a pre-existing list. 
        """
        dlist=string.split('\t') #tab delimited
        dlist=np.array([float(d) for d in dlist if len(d)>0])
        if len(dlist)>0:
            dataset.append(dlist)
        return len(dlist), dataset


    def import_data(self):
        """
        Extracts the data from the file and stores it in the class properties.
        """
        data = [] # initialize this as a list
        header = [] # initialize this as a list

        # go through the fle
        with open(self.f, 'rb') as f:
            for line in f:
                line = line.decode('latin-1')
                line = line.strip() # remove whitespace
                if line[0] == '#': # line is in the header
                    header.append(line)
                else: # line not in header
                    length, data = self.tsvtoarr(line, data)
                    # data is a list with the form [[wavenums],[point1_y,point1_x,[point1_intensities]],...,[pointN_y,pointN_x,[pointN_intensities]]]

        # sort out the components of the data list

        # check to see if you likely have eV or cm^-1
        if np.amin(np.asarray(data[0])) < 6:
            self.energy = np.asarray(data[0])
            self.wn = 1 / ((1240 / self.energy) * 1E-7)
        else:
            self.wn = np.asarray(data[0])
            self.energy = 1240 / ((1 / self.wn) * 1E7)
        self.position = np.zeros((len(data)-1), dtype=[('x', 'f4'), ('y', 'f4')])
        self.intensity = np.zeros((len(data)-1, len(self.wn)), dtype=np.float32)
     
        for i in range(1,len(data)):
            # check that this is a full map and not a line scan
            if len(data[i]) == (len(data[0]) + 2):
                self.position[i-1] = (data[i][1], data[i][0]) # x, y
                self.intensity[i-1] = data[i][2::]
            else:
                self.position[i-1] = (data[i][0], 0) # x, y
                self.intensity[i-1] = data[i][1::]



class RamanSpot():

    def __init__(self, file):
        self.f = file # filepath
        self.filename = self.get_filename()
        self.wn = None # ndarray of wavenumbers, shape (N, 1)
        self.energy = None # ndarray of energy, shape (N, 1)
        self.intensity = None # ndarray of peak intensities, shape (N, 1)
        self.import_data()


    def get_filename(self):
        """
        Takes in the entire path to a file. Returns just
        the filename as a string (no preceeding part of the path).
        """
        return (os.path.split(self.f)[1])[:-4]


    def import_data(self):
        """
        Extracts the data from the file. 
        """
        x_data = [] # initialize this as a list
        int_data = [] # initialize this as a list
        header = [] # initialize this as a list

        # go through the fle
        with open(self.f, 'rb') as f:
            for line in f:
                line = line.decode('latin-1')
                line = line.strip() # remove whitespace
                if line[0] == '#': # line is in the header
                    header.append(line)
                else: # line not in header
                    parts = line.split(' ') # split by space
                    x_data.append(parts[0])
                    int_data.append(parts[1])

        # turn the components into arrays
        # check to see if likely eV or cm^-1
        if np.amin(np.asarray(x_data)) < 5:
            self.energy = np.asarray(x_data)
            self.wn = 1 / ((1240 / self.energy) * 1E-7)
        else:
            self.wn = np.asarray(x_data)
            self.energy = 1240 / ((1 / self.wn) * 1E7)

        self.intensity = np.asarray(int_data)

     
