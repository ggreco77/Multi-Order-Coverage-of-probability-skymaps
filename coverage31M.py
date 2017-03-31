# -*- coding: utf-8 -*-

try:
   import cPickle as pickle
except:
   import pickle

import fileinput
import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import matplotlib.dates as mdates

import dateutil

from Tkinter import *
import tkMessageBox
import tkFont

from math import cos, sin, acos, asin, atan, degrees, radians

import astropy
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import TimeDelta
#from astroquery.vizier import Vizier
from astropy.table import Table

import numpy as np
import healpy as hp # fatto i moc puoi togliere
#from scipy.stats import norm

from mocpy import MOC

from aladinSAMP import AladinViaSAMP, AladinScriptCommands 
samp = AladinViaSAMP() # potresti toglierlo, ereditato da aladin cone jupyter
aladin = AladinScriptCommands()

from user_values import UserValues
from lvc_skymap import LVCskymap
from query import Query
from airmass import Airmass
from moon import Moon

from utils import Utils
Utils.create_folders(folders=["Queries", "Coords", "FoV"])
Utils.load_user_fov("GWsky_fov.vot")


class ShowSkyCoverage(object): 
    """Moving the FoV-footprint coverage by choosing a starting pointing."""
    
    SHIFT_CORRECTION = 0.00001  # A shift correction of 0.00001 is added
                                 # --> to escape math error during the FoV sequence

    def __init__(self, infile_coords='GWsky_coords'):
        """Creating a class in which the instance attributes are based on the dictionary
       "GWsky_coords" by default. GWsky_coords is created and pickled by the "user_values"
       module and will be deleted when the "ShowSkyCoverageGUI" will be closed.
       It contains the keys: "ra", "dec". ra and dec represents the central location of a FoV. 
        
        Starting sky coordinates:
             self.input_ra: right ascension [deg]
             self.input_dec: declination [deg]
         """
       
        self.infile_coords = infile_coords
        self.entries_GWsky_new =[] # new entries during the FoV sequence
        
        self.user = UserValues() # compositions
        self.lvc = LVCskymap()
        self.query = Query()
        self.airmass = Airmass()
        self.moon = Moon()
        
        with open(infile_coords, 'rb') as data:  
            coords_GWsky = pickle.load(data)
            
        for k, v in coords_GWsky.items():          
            setattr(self, k, v)
            
            
    def ra0ra1(self, A, dec0, dec1):
        """From the angular distance:
           cos(A) = sin(Dec1)sin(Dec2)+cos(Dec1)cos(Dec2)cos(ra1-ra2) --> 
           cos(ra1-ra2) = [cos(A)-sin(dec0)sin(dec1)]/[cos(dec0)cos(dec1)]."""

        dec0, dec1, A = radians(dec0),  radians(dec1), radians(A)
        cos_ra0_ra1 = ( cos(A)-sin(dec0)*sin(dec1) )/( cos(dec0)*cos(dec1) )
        ra0ra1 = degrees( acos(cos_ra0_ra1) )

        return  round(ra0ra1, 5)         
    
    def __updating_center_coords(self, ra, dec):
        """Getting/Updating FoV-center (ra, dec) in the dict named by default "GWsky_coords".
.          For each tile across the sky the file is updated."""""

        with open('GWsky_coords', 'rb') as data:
            coords_GWsky = pickle.load(data)
            
        coords_GWsky['input_ra'], coords_GWsky ['input_dec'] = ra, dec

        with open('GWsky_coords', 'wb') as data:
            pickle.dump(coords_GWsky, data)
   
    def __are_all_same(self, items):
        """Check if all elements of a list are the same."""
        
        return all(x == items[0] for x in items)

    def __fov_stats(self, ra, dec, table, integrated_prob, distance_grid, ansatz_distribution, moon_illumination):    
        """Managing the descriptive statistic window. If the airmass is outside the default range [airmass_min, airmass_max]
            the window is not opened otherwise the window is shown."""
        
        airmass_values, datestrings = self.airmass.airmass_step(ra, dec) 
        same = self.__are_all_same(airmass_values) 
                                                 
        if same == True:
            tkMessageBox.showerror('Warning',"airmass outside the range of {1 - 5.8}")
            aladin.remove_FoV(ra, dec) 
            aladin.remove("C_" + str( ra ) + "/" + str( dec ))  
        else:
            time_step = [dateutil.parser.parse(s) for s in datestrings]
            self.__updating_center_coords(ra,dec)  
            fov_statistics = FoVstatistics()                 
            fov_statistics.plot_stats(time_step, airmass_values,
                                      ra, dec,
                                      table,
                                      integrated_prob,
                                      distance_grid, ansatz_distribution, moon_illumination)
        print airmass_values, datestrings
    
    def update_pointings_file(self, infile, ra, dec, prob_fov):
         """The central location (ra[deg], dec[deg]) and the integrated probability of
             a selected FoV are saved locally in an external file.
             By default the file is named "GWsky_pointings.txt"."""
           
         with open(infile, 'a') as pointing:
             pointing.write(str(ra) + ' ' + str(dec)+ ' '+ str(prob_fov)+'\n')

    def __query_shape(self, ra, dec, fov_shape):
        """Return the catalog query according with the defined-user FoV shape:
                   (1) box and (2) circle. """
        
        if self.user.get_fov_shape() != 2:  # box FoV
                    query_result = self.query.query_box(
                       ra, dec, self.user.get_fov_width(), self.user.get_fov_height(), self.user.get_catalog(),
                       self.user.get_column_1(), self.user.get_column_2(),
                       self.user.get_filter_1(), self.user.get_filter_2())               
        else: # circle FoV
            query_result = self.query.query_circle(
               ra, dec, self.user.get_fov_radius(), self.user.get_catalog(),
               self.user.get_column_1(), self.user.get_column_2(),
               self.user.get_filter_1(), self.user.get_filter_2())
            
        return query_result

    def __prob_shape(self, ra, dec, fov_shape):
        """Return the integrated probability according with the defined-user FoV shape:
                   (1) box and (2) circle."""
        
        if self.user.get_fov_shape() !=2:  # box FoV
            prob_fov = self.lvc.prob_in_box(
               ra, dec, self.user.get_fov_width(), self.user.get_fov_height())                
        else: # circle FoV
            prob_fov = self.lvc.prob_in_circle(
               ra, dec, self.user.get_fov_radius())
            
        return  prob_fov
                  
    def pick_coverage(self, ra, dec):
        """Setting GWsky: with statistic window (A); without statistic window (D)."""                                      

        if self.user.get_GWsky_basic() == "A":  # full version 

            query_result = self.__query_shape(ra, dec, self.user.get_fov_shape())
            prob_fov = self.__prob_shape(ra, dec, self.user.get_fov_shape())
            moon_illumination =  self.moon.illumination()              # moon_illumination
                                                                             # moon_dist
            # TEST   
            fov_sep = Utils.separation(self.input_ra, self.input_dec, 
                                        ra, dec)
            print ('The distance between 2 consecutive FoV centers is', fov_sep.round(6))
            
            r, dp_dr = self.lvc.conditional_distance_fov_center(ra, dec) 
            self.__fov_stats(ra, dec, query_result, prob_fov, r, dp_dr, moon_illumination) # add moon ill, moon dist
            
        elif self.user.get_GWsky_basic() == "D":  # basic version-> no Stats win
           
            prob_fov = self.__prob_shape(ra, dec, self.user.get_fov_shape())

            # TEST  
            #Utils.separation(self.input_ra, self.input_dec, ra, dec)                                       

            self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov)
            
    def intercardinal_distance(self, ra, dec, shift_up_down, shift_right_left):
        """Moving from the fixed cardinal direction using the bi-directional windows;
           shift_up_down ↕ and/or shift_right_left ↔."""

        if shift_right_left > 0:
           shift_east_west = self.ra0ra1((shift_right_left-self.SHIFT_CORRECTION),
                                                   (dec + self.user.get_fov_height() + shift_up_down),
                                                   (dec + self.user.get_fov_height() + shift_up_down))
           dist = ra + shift_east_west 
         
        elif shift_right_left < 0 :
           shift_east_west = self.ra0ra1((shift_right_left + self.SHIFT_CORRECTION),
                                                   (dec + self.user.get_fov_height() + shift_up_down),
                                                   (dec + self.user.get_fov_height() + shift_up_down))
           dist = ra - shift_east_west
         
        else:
           dist = ra

        return dist 

    def load_entries(self, infile_entries):
        """Opening the file in which the input entries are stored: ra_1 dec_1 ra_2 dec_2 ... ra_n dec_n
           By default the file is named "GWsky_entries". "GWsky_entries" is created from the
           "show_starting_fov" method in "StartingFoV" class.
           A error message invites users to press the "Start FoV" button before carrying out any"""
       
        try:
            with open(infile_entries, 'rb') as data:
                entries_GWsky = pickle.load(data)
                return entries_GWsky
        except IOError as io_error:
            message = "Press the button 'Start FoV' and choose a starting position; \
                        by default the sky position of the max probability pixel is shown"
            tkMessageBox.showerror ('Error', message)
                         
    def north(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in North direction from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔ """

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            north_pointing = [(dist),
                               (float(dec_start) + self.user.get_fov_height() + shift_up_down)]
             
            ra, dec = north_pointing[0], north_pointing[1]

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
            
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def south(self, shift_up_down=0, shift_right_left=0):    
        """Moving the FoV tiles in South direction from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔"""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            dist = self.intercardinal_distance(float(ra_start), float(dec_start),
                                                 shift_up_down, shift_right_left)
            south_pointing = [(dist), (float(dec_start) - self.user.get_fov_height() - shift_up_down)]
                    
            ra, dec = south_pointing[0], south_pointing[1]
            
            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
                                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)    
          
    def east(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in East direction  from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔.
           A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):              
            ra_distance = self.ra0ra1((self.user.get_fov_width() - self.SHIFT_CORRECTION + shift_right_left),
                                        float(dec_start), float(dec_start))
                
            east_pointing = [(float(ra_start) + ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = east_pointing[0], east_pointing[1]

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)           

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
                   
    def west(self, shift_up_down=0, shift_right_left=0):
        """Moving the FoV tiles in West direction  from input sky-position(s).
           The related bidirectional button permits small shifts from such pre-defined direction:
           shift_up_down ↕ and/or shift_right_left ↔.
            A shift correction of 0.00001 is added to escape math error."""
        
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]  
      
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
               
            ra_distance = self.ra0ra1((self.user.get_fov_width() - self.SHIFT_CORRECTION + shift_right_left),
                                      float(dec_start), float(dec_start))

            west_pointing = [(float(ra_start) - ra_distance), (float(dec_start) + shift_up_down)]
            ra, dec = west_pointing[0], west_pointing[1] 

            aladin.get_FoV(ra, dec)
            self.pick_coverage(ra, dec)
            
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
        

class ShowSkyCoverageGUI(Toplevel):
    """
        The GUI consists of 9 buttons; 4  directional buttons to move the FoVs
        in cardinal directions (North, South, East, West) 4  buttons to shift the FoVs
        from a consecutive cardinal direction (↕, ↕, ↔, ↔) and 1 button to get a
        new FoV position (Start FoV). ***Input values in deg***.
    """
    
    def __init__(self, tkMainWin):
        frame = Frame(tkMainWin, border=9, bg="dim grey")
        frame.pack()     

        self.B02 = Button(frame,text="N", )   
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.grid(row=0, column=2)
        
        self.B12 = Button(frame,text="↕↔",
                          fg="grey")  
        self.B12.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B12.grid(row=1, column=2)
  
        self.B30 = Button(frame, text="E", )  
        self.B30.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B30.grid(row=3,column=0)

        self.B31 = Button(frame,text="↕↔",
                          fg="grey")   
        self.B31.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B31.grid(row=3, column=1)

        self.B32 = Button(frame, text="Start FoV",
                          fg="red4",)
        self.B32.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B32.grid(row=3, column=2)

        self.B33 = Button(frame,text="↕↔",
                          fg="grey")  
        self.B33.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B33.grid(row=3, column=3)

        self.B34 = Button(frame, text="W",) 
        self.B34.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B34.grid(row=3, column=4)

        self.B42 = Button(frame,text="↕↔",
                          fg ="grey") 
        self.B42.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B42.grid(row=4, column=2)
        
        self.B52 = Button(frame,text="S", ) 
        self.B52.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B52.grid(row=5, column=2)

                                 ## Adjustments Btns
        self.B60 = Button(frame,text="↞",
                          fg ="grey",pady=1,) 
        self.B60.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B60.grid(row=6, column=0)

        self.B61 = Button(frame,text="↠",
                          fg ="grey",pady=1,) 
        self.B61.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B61.grid(row=6, column=1)

        self.B62 = Button(frame,text="✓ Accept",fg ="green4",
                          pady=1) 
        self.B62.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B62.grid(row=6, column=2)

        self.B63 = Button(frame,text="↟",
                          fg ="grey",pady=1) 
        self.B63.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B63.grid(row=6, column=3)

        self.B64 = Button(frame,text="↡",
                          fg ="grey",pady=1) 
        self.B64.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B64.grid(row=6, column=4)

                                   ## Folder
        self.B72 = Button(frame,text="⏩▶ Folder",
                          fg ="gold4",pady=1) 
        self.B72.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B72.grid(row=7, column=2)
        

    # Actions    
    def clicked(self, event):
        """Moving the user-defined FoV footprint."""
        
        run_sequence = ShowSkyCoverage()
        
        if event.widget == self.B02:
            run_sequence.north()         # north

        if event.widget == self.B12:
            move_fov = ShiftFoV()
            move_fov.north_shift()       # ↕↔
                      
        if event.widget == self.B30:
            run_sequence.east()          # east

        if event.widget == self.B31:
            move_fov = ShiftFoV()
            move_fov.east_shift()        # ↕↔
            
        if event.widget == self.B32: 
            new_starting_fov = StartingFoV()  # start FoV
            new_starting_fov

        if event.widget == self.B33:
            move_fov = ShiftFoV()
            move_fov.west_shift()        # ↕↔
            
        if event.widget == self.B34:   
            run_sequence.west()          # west
            
        if event.widget == self.B42:
            move_fov = ShiftFoV()
            move_fov.south_shift()       # ↕↔
            
        if event.widget == self.B52:    
            run_sequence.south()         # south

        if event.widget == self.B60:     # ↞        
            adj = Adjustments()
            adj.adj_east()
            
        if event.widget == self.B61:     # ↠       
            adj = Adjustments()
            adj.adj_west()            
                
        if event.widget == self.B63:     # ↟          
            adj = Adjustments()
            adj.adj_north()

        if event.widget == self.B64:     # ↡           
            adj = Adjustments()
            adj.adj_south()         
               
        if event.widget == self.B62:     # ✓ Accept
            adj = Adjustments()
            adj.adj_accept()

        if event.widget == self.B72:     # ▶ Folder
            Utils.move_to_folder(planes=['Q:*','P:*'],
                                 folders=['Queries','FoV'])
            

class Adjustments(ShowSkyCoverage):
    """Adjustments FoV position."""

    def __init__ (self):

       ShowSkyCoverage.__init__(self)
       
       self.shift_up = 0.15     # default adjustments   (up)
       self.shift_down = 0.15   #      "               (down)
       self.shift_left = 0.15   #      "               (left)
       self.shift_right = 0.15  #      "               (right)
       
    def adj_north(self):
        """Adjustments FoV position -> north direction"""
            
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
            
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
               
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            dist = self.intercardinal_distance(ra_start, dec_start,
                                               self.shift_up, shift_right_left=0)
            north_adj = [(dist),
                         (dec_start + 0 + self.shift_up)]
             
            ra, dec = north_adj[0], north_adj[1]
                
            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))
                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)
            
            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))  

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_south(self):
        """Adjustments FoV position -> south direction"""
         
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
            
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
               
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            dist = self.intercardinal_distance(ra_start, dec_start,
                                               self.shift_down, shift_right_left=0)
            south_adj = [(dist),
                         (dec_start + 0 - self.shift_down)]
             
            ra, dec = south_adj[0], south_adj[1]
                
            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))
                
            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)
            
            #aladin.remove("Q:"+str(ra_start)+"/"+str(dec_start))
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))  

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_east(self):
        """Adjustments FoV position -> east direction"""

        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)

            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
            
            ra_distance = self.ra0ra1((0 - self.SHIFT_CORRECTION + self.shift_left),
                                        float(dec_start), float(dec_start))
                          
            east_adj = [(float(ra_start) + ra_distance), (float(dec_start) + 0)]
            ra, dec = east_adj[0], east_adj[1]

            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))       

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))  

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)
            
    def adj_west(self):
        """Adjustments FoV position -> weast direction"""
         
        entries_GWsky = self.load_entries("GWsky_entries")
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]

        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)

            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
             
            ra_distance = self.ra0ra1((0 - self.SHIFT_CORRECTION + self.shift_right),
                                        float(dec_start), float(dec_start))
            
            aladin.select("P:"+str(ra_start) + ',' + str(dec_start))
                
            west_adj = [(float(ra_start) - ra_distance), (float(dec_start) + 0)]
            ra, dec = west_adj[0], west_adj[1]

            aladin.set_target(ra, dec)
            aladin.set_plane_id("P:"+str(ra) + ',' + str(dec))       

            new_sky_pos = [ra,dec] # cycle variables
            self.entries_GWsky_new.extend(new_sky_pos)

            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                  ra=str(ra_start), dec=str(dec_start))  

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

    def adj_accept(self):
        """Confirming the adjustments FoV position -> open statistic win."""
         
        entries_GWsky = self.load_entries("GWsky_entries")        
        fov_center_ra, fov_center_dec = entries_GWsky[0::2], entries_GWsky[1::2]
                   
        for ra_start, dec_start in zip (fov_center_ra, fov_center_dec):
            ra_start, dec_start = float(ra_start), float(dec_start)
                    
            self.pick_coverage(float(ra_start), float(dec_start))

            new_sky_pos = [ra_start,dec_start]
            self.entries_GWsky_new.extend(new_sky_pos)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky_new, data)

            
class StartingFoV(Toplevel):
    """Starting a sequence from a list of FoV(s). The window contains 1 keyboard entries and 3 Buttons.

        entry:
             ra_1 dec_1 ra_2 dec_2 ra_3 dec_3 ... ra_n dec_n [deg]
             By default: sky coords of maximum probability pixel

        Btns:
           Show : draw user-defined FoV footprint(s) in Aladin Plane(s)
           No show : no draw user-defined FoV footprint(s) in Aladin Plane(s)
           Close : close the widget
        """
    
    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()

        # putting the entry value(s) in a list
        self.entries_GWsky=[]  

        self.wait_visibility()
        self.wm_attributes('-alpha', 0.8) # transparency

        self.title(" Starting FoV")
        self.attributes("-topmost", True)
        
        self.label_1 = Label(self, text="RA (°) DEC (°)", bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, pady=0)

        # default: sky coords of maximum probability pixel
        fov_coords = str(self.user.get_ra_max_pixel()), str(self.user.get_dec_max_pixel()) 
        
        max_pixel_default = StringVar(self, value=fov_coords) 
        self.entry_1 = Entry(self, width=30, justify=CENTER,
                             textvariable=max_pixel_default)

        self.entry_1.grid(row=0, padx=15, column=1)

        self.entryScroll = Scrollbar(self, orient=HORIZONTAL,
                                     command=self.__scrollHandler)
        self.entryScroll.grid(row=1, column=1, sticky=E+W)
        self.entry_1['xscrollcommand'] = self.entryScroll.set

        #Btn

        self.show = Button(self, text='Show',
                           command=self.show_starting_fov)
        self.show.grid(column=2, row=0, sticky=W, padx=2, pady=5)
        
        self.checkbox = Button(self, text="Not show",      
                               command=self.no_show_starting_fov)
        self.checkbox.grid(column=3,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Obs",  fg='dark green',
                            command=self.obs)  
        self.close.grid(column=4,row=0, sticky=W, padx=2, pady=5) 
        
        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=5,row=0, sticky=E, padx=2, pady=5)

    #Actions
    def __scrollHandler(self, *L):
        """Scroll entry in starting FoV window."""
       
        op, howMany = L[0], L[1]

        if op == 'scroll':
            units = L[2]
            self.entry_1.xview_scroll(howMany, units)
        elif op == 'moveto':
            self.entry_1.xview_moveto(howMany)
            
    def __split_entries(self):
        """Splitting the entries in ra and dec; # odd: ra and # even: dec."""
        
        current_fov_coords = self.entry_1.get().replace(';',' ').replace(',',' ').split()
        fov_center_ra = current_fov_coords[0::2]
        fov_center_dec = current_fov_coords[1::2]

        return current_fov_coords, fov_center_ra, fov_center_dec
    
    def show_starting_fov(self):
        """Drawing the FoV footprint(s) in the Aladin plane(s).
         By default: sky coords (ra[deg], dec[deg]) of maximum probability pixel."""   
   
        show_sky_coverage = ShowSkyCoverage()
        
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()
        
        try:
            for ra_starting, dec_starting in zip (fov_center_ra, fov_center_dec):
                aladin.get_FoV(float(ra_starting), float(dec_starting))
                show_sky_coverage.pick_coverage(float(ra_starting), float(dec_starting))
                              
        except ValueError as value_error:
            tkMessageBox.showerror ('Error', value_error)
          
        self.entries_GWsky.extend(current_fov_coords)

        with open('GWsky_entries', 'wb') as data:
            pickle.dump(self.entries_GWsky, data)
            
        self.entries_GWsky=[] # re-init.
        
    def no_show_starting_fov(self):
        """No Draw the FoV footprint(s) in the Aladin plane(s);
               useful to re-initialize the sequence."""

        self.entries_GWsky=[] # re-init.     
        current_fov_coords, fov_center_ra, fov_center_dec = self.__split_entries()

        self.entries_GWsky.extend(current_fov_coords)

        with open('GWsky_entries', 'wb') as data:
            return pickle.dump(self.entries_GWsky, data)

    def obs(self):
        """Inizializing observability window."""
        
        observability = Observability()
     
    def close_window(self):
        self.destroy()


class Observability(Toplevel):
    """GUI to create a MOC region in which the airmass is less than a value defined by users."""

    def __init__(self):
        Toplevel.__init__(self, border=8, bg="slate grey")
        
        self.user = UserValues()

        self.wait_visibility()
        self.wm_attributes('-alpha',0.8) # transparency

        self.title("Observability" + "@" + self.user.get_obs_time())
        self.attributes("-topmost", True)
        
        # first label
        self.label_1 = Label(self, text="Show the region in the",
                             bg="slate grey")
        self.label_1.grid(row=0, column=0, sticky=E, padx=0)

        moc_value = 90  # default     
        moc_default = StringVar(self, value=moc_value)
        
        self.entry_percentage = Entry(self, width=5, justify=CENTER,
                             textvariable=moc_default)
        self.entry_percentage.grid(row=0, padx=2, column=1)

        # second label
        self.label_2 = Label(self, text="% MOC in which the airmass is ≤",
                             bg="slate grey")
        self.label_2.grid(row=0, column=2, sticky=E, pady=0)

        airmass_value = "2.5" # default
        airmass_default = StringVar(self, value=airmass_value)
        
        self.entry_airmass = Entry(self, width=5, justify=CENTER,
                             textvariable=airmass_default)
        self.entry_airmass.grid(row=0, padx=2, column=3)

        #Btn
        self.show = Button(self, text='Show',
                           command=self.moc_obs)
        self.show.grid(column=4, row=0, sticky=W, padx=2, pady=5)
        
        self.moon = Button(self, text="Moon",
                            command=self.close_window)  
        self.moon.grid(column=5,row=0, sticky=W, padx=2, pady=5) 
        
        #self.forward = Button(self, text=">>",      
        #                       command=self.close_window)
        #self.forward.grid(column=6,row=0, sticky=E, padx=2, pady=5)

        self.close = Button(self, text="Close",
                            command=self.close_window)  
        self.close.grid(column=7,row=0, sticky=E, padx=2, pady=5)


    #Actions
    def moc_obs(self):
        """Creating a MOC region in which the airmass is less than a value defined by users."""

        percentage = float(self.entry_percentage.get())/100.0

        hpx = hp.read_map(self.user.get_skymap(), verbose = False)
        
        sort = sorted(hpx, reverse = True)
        cumsum = np.cumsum(sort)
        index, value = min(enumerate(cumsum), key = lambda x: abs(x[1] - percentage))
        value = round(value, 1) # value --> threshold
        print percentage
        print index, value

        # finding ipix indices confined in a given percentage 
        index_hpx = range(0, len(hpx))
        hpx_index = np.c_[hpx, index_hpx]

        sort_2array = sorted(hpx_index, key = lambda x: x[0], reverse = True)
        value_contour = sort_2array[0:index]

        j = 1 
        table_ipix_contour = [ ]

        for i in range (0, len(value_contour)):
            ipix_contour = int(value_contour[i][j])
            table_ipix_contour.append(ipix_contour)
          
     
        # from index to polar coordinates
        theta, phi = hp.pix2ang(self.user.get_nside(), table_ipix_contour)

        # converting these to right ascension and declination in degrees
        ra = np.rad2deg(phi)
        dec = np.rad2deg(0.5 * np.pi - theta)

        #print ra

        # airmass calculation
        obs_time = Time(self.user.get_obs_time())

        observatory = astropy.coordinates.EarthLocation(
           lat=self.user.get_latitude()*u.deg, lon=self.user.get_longitude()*u.deg,height=self.user.get_altitude()*u.m)

        sky_coord = SkyCoord(ra = ra*u.deg,dec=dec*u.deg, frame='icrs')

        altaz = sky_coord.transform_to(AltAz(obstime=self.user.get_obs_time(),
                                             location=observatory))

        airmass_values = altaz.secz
        print airmass_values
        # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
        contour_ipix = Table([ ra, dec, airmass_values ], names = (
           'RA[deg]', 'DEC[deg]', 'airmass'), meta = {'ipix': 'ipix table'})

        mask = (contour_ipix['airmass']) >1 # clear

        obs1=contour_ipix[mask]

        mask2 = (obs1['airmass']) <= float(self.entry_airmass.get())  # user values

        obs = obs1[mask2]

        # setting MOC order
        from math import log
        moc_order = int(log( self.user.get_nside(), 2))

        # creating a MOC map from the contour_ipix table
        moc = MOC.from_table( obs, 'RA[deg]', 'DEC[deg]', moc_order )

        # writing MOC file in fits
        moc.write( 'obs_airmass_'+self.entry_airmass.get()+'MOC_'+str(percentage), format = 'fits' )
        print type(self.entry_airmass.get())

    def close_window(self):
        self.destroy()

        
class ShiftFoV(Toplevel):
    """Shifting the FoV footprint(s) from a consecutive cardinal direction (↕, ↕, ↔, ↔);
       The widget contains 2 entries and 2 Buttons."""
    
    def __init__(self):
        Toplevel.__init__(self, border=7, bg="slate grey")
        self.attributes("-topmost", True)
        self.wait_visibility()
        self.wm_attributes('-alpha', 0.8)     
        
        self.label_3 = Label(self, text="↕ (°)",bg="slate grey")  
        self.entry_3 = Entry(self, width=6, justify=CENTER)
        self.label_3.grid(row=0, sticky=E) 
        self.entry_3.grid(row=0, column=1)

        self.label_4 = Label(self, text="↔ (°)",bg="slate grey")  
        self.entry_4 = Entry(self, width=6, justify=CENTER)
        self.label_4.grid(row=0,column=3) 
        self.entry_4.grid(row=0, column=4)

        self.close = Button(self, text="Close",
                            command = self.close_window)  
        self.close.grid(column=4,row=2)
    
    def north_shift(self):
        self.title(" Shifting - North")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_north)
        self.checkbox.grid(column=3,row=2)

    def south_shift(self):
        self.title(" Shifting - South")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_south)
        self.checkbox.grid(column=3,row=2)

    def east_shift(self):
        self.title(" Shifting - East")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_east)
        self.checkbox.grid(column=3,row=2)

    def west_shift(self):
        self.title(" Shifting - West")
        self.checkbox = Button(self, text="OK",
                               command=self.shift_west)
        self.checkbox.grid(column=3,row=2)
        
    # Actions
    def shift_north(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.north(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_south(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.south(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_east(self):             
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.east(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def shift_west(self):              
        try:
            shift_up_down = float(self.entry_3.get())
            shift_right_left = float(self.entry_4.get())
            shift = ShowSkyCoverage()
            shift.west(shift_up_down, shift_right_left)
        except ValueError as value_error:
            tkMessageBox.showerror ('Error',value_error)

    def close_window(self):
        self.destroy()

class FoVstatistics(Toplevel, ShowSkyCoverage):
    """ FoV statistics window consists of 3 buttons:
         (1) Confirm the FoV-footprint
         (2) Delete FoV-footprint
         (3) Zoom in FoV. """
    
    def __init__(self):
        Toplevel.__init__(self)
        ShowSkyCoverage.__init__(self)

        self.title("FoV statistics")
        self.attributes("-topmost", True)
              
        self.B00 = Button(self,text="Confirm the pointing in the GWsky_pointings txt file")  
        self.B00.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B00.pack(side=TOP,fill=BOTH)

        self.B01 = Button(self,text="Delete the FoV")  
        self.B01.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B01.pack(side=TOP,fill=BOTH)

        self.B02 = Button(self,text="Zoom in the FoV")  
        self.B02.bind("<ButtonRelease-1>",
                      self.clicked)
        self.B02.pack(side=TOP,fill=BOTH)

        # default
        fov_center_ra_dec = str(self.input_ra), str(self.input_dec), 'on fly notes, NO DELETE OR CHANGE THE COORDS!'       
        current_fov = StringVar(self, value=fov_center_ra_dec)  
        self.entry_current_fov = Entry(self, width=30, justify=LEFT,
                                       textvariable=current_fov)
        
        #self.entry_current_fov.bind("<Key>", lambda e: "break") # saving entries
        self.entry_current_fov.pack(side=TOP,fill=BOTH)

    def __rm_from_stack(self, ra, dec):
        """Removing from Aladin stack the associated planes."""
        
        aladin.remove_FoV(ra, dec)       
        aladin.remove("Q:"+ ra +"/"+ dec)           
        #aladin.remove("C_" + ra+ "/" + dec)

    def clicked(self, event):
        """Retain or delate the FoV-footprint(s).
            The FoV-center positions are saved in "Pointings txt" """      
        
        if event.widget == self.B00: # Retain and Close the FoV.            
            self.destroy()
            
        if event.widget == self.B01: # Delete FoV
            current_fov_coords = self.entry_current_fov.get().split() # getting entries
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]
            
            self.__rm_from_stack(current_fov_ra,current_fov_dec)
            
            Utils.delete_pointing(infile="GWsky_pointings.txt",
                                 ra=str(current_fov_ra), dec=str(current_fov_dec))           
            self.destroy()
       
        if event.widget == self.B02:  # Zoom in the  FoV
            aladin.location(str(self.input_ra), str(self.input_dec))             
            aladin.zoom('1x')
            
    # Plots        
    def plot_stats(self, time_step, airmass_values, ra, dec, table_q,
                   prob_fov, r, dp_dr, moon_illumination):
        """Showing the plots in the FoV statistic window."""
        
        f = Figure(figsize=(9, 5.2), facecolor='white')
        f.subplots_adjust(left=.13, bottom=.16, right=.93, top=.84, wspace=.26, hspace=.3)

        def airmass_subplot():
            """SubPlot Airmass."""
            
            airmass_plt = f.add_subplot(223)   
            ax=plt.gca()
            ax.xaxis.set_major_formatter(
                mdates.DateFormatter('%H:%M:%S'))

            airmass_plt.set_ylabel('airmass')
            airmass_plt.set_xlabel('Universal Time')
            airmass_plt.invert_yaxis()
            airmass_plt.grid(True)

            airmass_plt.plot(time_step, airmass_values, "o-", linestyle="--")
            f.autofmt_xdate()

        airmass_subplot = airmass_subplot()

        def subplot_cat_distribution():
            """SubPlot catalogue distributions."""

            # split entry: ra and dec
            current_fov_coords = self.entry_current_fov.get().split() 
            current_fov_ra, current_fov_dec = current_fov_coords[0], current_fov_coords[1]

            try:
                for table_name in table_q.keys():
                    table = table_q[table_name]
                
                table.write("GWsky_query_items", format = 'votable', overwrite = True)
                
                query_catalog = f.add_subplot(222)

                # column 1
                query_catalog.hist(table[self.user.get_column_1()].quantity)         
                query_catalog.set_xlabel('cat: ' + self.user.get_catalog() + ' ' +
                                         'col: ' + self.user.get_column_1())
                query_catalog.set_ylabel('Count')                     

                # column 2
                query_catalog = f.add_subplot(224) 
                query_catalog.hist(table[self.user.get_column_2()].quantity)                   
                query_catalog.set_xlabel('cat: ' + self.user.get_catalog() + ' ' +
                                         'col: ' + self.user.get_column_2())
                query_catalog.set_ylabel('Count')                                
                
            except UnboundLocalError as unbound_local_error:
                #print ('No Galaxies in the selected FoV')
                tkMessageBox.showerror(' ', unbound_local_error)
            except KeyError as key_error:
                #print ('No catalog column:', key_error)
                tkMessageBox.showerror(' ', key_error)
            finally:
                pass
            
        subplot_cat_distribution = subplot_cat_distribution()       

        def subplot_cond_distance():
            """Conditional Distance Distribution Along a Line of Sight (FoV center position)."""
            
            conditional_distance_line_sight = f.add_subplot(221) 
            conditional_distance_line_sight.plot(r, dp_dr)
            
            title_string = 'Conditional distance distribution \n along the FoV center'
            conditional_distance_line_sight.set_title(title_string,fontsize=10)
            conditional_distance_line_sight.set_xlabel('distance (Mpc)')
            conditional_distance_line_sight.set_ylabel('prob Mpc$^{-1}$')

        subplot_cond_distance = subplot_cond_distance()

        def draw_area():
            """Drawing Area of the FoV Statistic window."""

            fov_information_title = "FoV center (ra "+str(ra) + "  " +"dec "+ str(dec)+")" + "; " + "prob: " + str(prob_fov)+ \
                                    ";Moon(illumi.: " + str(moon_illumination) + "and dist.:)"    
            f.suptitle(fov_information_title, fontsize=10)
           
            canvas = FigureCanvasTkAgg(f, self) 
            canvas.show()
            canvas.get_tk_widget().pack(side=BOTTOM, fill=BOTH, expand=True)

            toolbar = NavigationToolbar2TkAgg(canvas, self)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=True)
            
        draw_area = draw_area()

        self.update_pointings_file("GWsky_pointings.txt", ra, dec, prob_fov)

        samp.send_file("GWsky_query_items")
        aladin.rename("Q:"+str(ra)+"/"+str(dec))

def on_closing():
    """Asking the closure of the coverage window. If "Quit" the files in the list "temp_files" are deleted.
          ***Improving with tempfile module***"""
    
    if tkMessageBox.askokcancel("Quit", "Do you want to quit?"):
        try:
            temp_files=["GWsky_entries", "GWsky_query_items", "GWsky_fov.vot",
                        "GWsky_config", "GWsky_coords"]
            for temp_file in temp_files:
               os.remove(temp_file)
        except OSError:
            pass
        mainWin.destroy()

        
# running
mainWin = Tk()

sscGUI = ShowSkyCoverageGUI(mainWin)
mainWin.title('GWsky')
mainWin.attributes("-topmost", True)

mainWin.wait_visibility(mainWin)
mainWin.wm_attributes('-alpha', 0.8) # transparency

mainWin.protocol("WM_DELETE_WINDOW", on_closing)
mainWin.mainloop()
