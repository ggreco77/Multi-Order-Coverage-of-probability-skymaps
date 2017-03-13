from __future__ import print_function

# Sending Aladin script commands

def send_script( script ):
    
    """
    
    Sending script commands to Aladin via SAMP
    
    """
    
    from astropy.vo.samp import SAMPIntegratedClient
     
    client = SAMPIntegratedClient()
    client.connect()

    params = {}
    message = {} 
    message[ "samp.mtype" ] = "script.aladin.send"
    message[ "samp.params" ] = { "script" : script }  

    client.notify_all( message )

    client.disconnect()

    
def cview( url ): 
    
    """
    
    Creation of view: url
    
    """
    cview_url = 'cview' + ' ' + url
    send_script( cview_url )

    
def rename ( plane ):
    
    """
    
    Rename plane

    """
    
    rename_plane = 'rename' + ' ' + plane
    send_script( rename_plane )
        
        
def get_hips ( catalog ): 
    
    """
    
    Call a remote image or tabular data server
    
    """

    hips = 'get' + ' ' + 'hips(' + catalog + ')'
    send_script ( hips )
    

def cview_plane ( plane ):
    
    """
    
    Creation of view: plane
      
    """
    
    show_plane = 'cview'  + ' "'+plane+'"'
    send_script ( show_plane )
    
    
def hide ( plane ):
    
    """
    
    Hide plane
     
    """
    
    hide_plane = 'hide'  + ' '+ plane
    send_script ( hide_plane )
    
    
def draw_line ( line_values ):
    
    """
    
    Graphical overlay commands: drawing a line
    
    """
    
    draw_line = 'draw' + ' ' + 'line' + ' ' + line_values
    send_script( draw_line )


def draw_newtool (name):
    
    """
    
    Graphical overlay commands:
    creating  a new drawing plane

    """
    
    draw_newtool = 'draw' + ' ' + 'newtool' + ' ' + name
    send_script( draw_newtool )
    
    
def rm_all():
    
    """
    
    Removing all planes
    
    """
    
    rm_all = 'rm -all'
    send_script( rm_all )
    
    
def send_file( infile ):
    
    """
    
    Sending file/table to Aladin plane via SAMP
        
    """

    global params
     
    from astropy.vo.samp import SAMPIntegratedClient
     
    client = SAMPIntegratedClient()
    client.connect()
    params = {}
    
    import sys
    import os.path
    
    if sys.version > '3':
        import urllib.parse
        params[ "url" ] = urllib.parse.urljoin( 'file:', os.path.abspath( infile ) )
    else:
        import urlparse
        params[ "url" ] = urlparse.urljoin( 'file:', os.path.abspath( infile ) )

    message = {}
    message[ "samp.mtype" ] = "image.load.fits"
    message[ "samp.params" ] = params

    client.notify_all( message )
    client.disconnect()
    
    
    ###############################
    def get_json_link( json_link ):

    """
    
    Plotting contour lines from a specific url
    
    """
    
    import json
    import sys
    import os.path
    
    # download the json file from "url" and save it locally under "contour.json"
    if sys.version < '3':
        import urllib
        jsonfile = urllib.URLopener()
        jsonfile.retrieve( json_link, "contour.json" )
    
    else:
        import urllib.request
        urllib.request.urlretrieve( json_link, "contour.json" )
     
    with open( 'contour.json' ) as data_file:
       data = json.load( data_file )

    contour_pieces = len( data[ 'contours' ] )

    percentile = ('10-percentile','20-percentile','30-percentile','40-percentile',
                  '50-percentile','60-percentile','70-percentile',
                  '80-percentile','90-percentile')

    for percentile_json in percentile:
       draw_newtool ( percentile_json )
       plot_contours_from_json( data, contour_pieces, percentile_json )
        
        
def plot_contours_from_json( data, contour_pieces, percentile_json ):

    """
    
    Managing the contour lines in a LVC json file
        
    """
    
    i = 0
    for i in range( 0, contour_pieces ):
        contour = data[ 'contours' ][ i ]
        percentile = contour[ 'name' ]

        if percentile == percentile_json:
            values = contour[ 'coords' ]

            # sending Aladin plane
            line = ( str( values ).replace('[' , '' ).replace(']' , '') )
            draw_line ( line )


def MOC_confidence_region( infile, percentage, short_name = ' ' ):
      
    """
    
    Multi-Order coverage map (MOC) of sky area enclosed within a contour plot
    at a given confidence level.
    
    Input:
         infile: healpix format
                 LVC probability sky map
         percentage: float
                  probability percentage of the enclosed area  
         short_name: str
                 output file name
     
    Output: fits format
                 MOC map named "short_name"_"percentage" 
                 
                 Remark: for json format change the statement
                 "moc.write(short_name+'_MOC_'+str(percentage), format='fits' )" -->  
                 "moc.write(short_name+'_MOC_'+str(percentage), format='json' )"        

    """
 
    import healpy as hp
    import numpy as np
     
    #reading skymap
    hpx = hp.read_map( infile, verbose = False )
    npix = len( hpx )
    nside = hp.npix2nside( npix )
 
    sort = sorted( hpx, reverse = True )
    cumsum = np.cumsum( sort )
    index, value = min( enumerate( cumsum ), key = lambda x: abs( x[1] - percentage ) )

    # finding ipix indices confined in a given percentage 
    index_hpx = range( 0, len( hpx ) )
    hpx_index = np.c_[ hpx, index_hpx ]

    sort_2array = sorted( hpx_index, key = lambda x: x[0], reverse = True )
    value_contour = sort_2array[ 0:index ]

    j = 1 
    table_ipix_contour = [ ]

    for i in range ( 0, len( value_contour ) ):
        ipix_contour = int( value_contour[i][j] )
        table_ipix_contour.append( ipix_contour )
          
     
    # from index to polar coordinates
    theta, phi = hp.pix2ang( nside, table_ipix_contour )

    # converting these to right ascension and declination in degrees
    ra = np.rad2deg( phi )
    dec = np.rad2deg( 0.5 * np.pi - theta )


    # creating an astropy.table with RA[deg] and DEC[deg] ipix positions
    from astropy.table import Table
    contour_ipix = Table([ ra, dec ], names = ('RA[deg]', 'DEC[deg]'), 
                         meta = {'ipix': 'ipix table'})
     
    
    # setting MOC order
    from math import log
    moc_order = int( log( nside, 2 ) )

    # creating a MOC map from the contour_ipix table
    moc = MOC.from_table( contour_ipix, 'RA[deg]', 'DEC[deg]', moc_order )

    # writing MOC file in fits
    moc.write( short_name + '_MOC_' + str( percentage ), format = 'fits' )

    # sending to Aladin plane
    send_file( short_name + '_MOC_' + str( percentage ) )
    cview( url = str( params[ 'url' ]) )
    rename ( plane = short_name  + '_MOC_' + str( percentage ) )
    
    
    #2. Multiscale meshes of gravitational-wave sky maps using MOC
    # selecting an event id (2015);
# http://www.ligo.org/scientists/first2years/
event_id = '18951'

# bayestar sky map
skymap_pipeline = 'bayestar'

# setting enclosed probability percentage 
prob_percentage = 0.9

# loading the simulated CBC event id (2015)
from astropy.utils.data import download_file

url_id = 'http://www.ligo.org/scientists/first2years/2015/compare/'+event_id+'/'+skymap_pipeline+'.fits.gz'
pipeline_event = download_file( url_id, cache = True, timeout = 300 )

# sending to Aladin plane
send_file ( pipeline_event )
rename ( skymap_pipeline + event_id )

# plotting contours from a specific url
from mocpy import MOC
get_json_link( 'https://losc.ligo.org/s/skymapViewer/json/skymaps/F2Y/'+event_id+'.json' )

# MOC extraction: 
#        area enclosed within a specific contour plot at a given confidence level
MOC_confidence_region( infile = pipeline_event, percentage = prob_percentage, 
                      short_name = skymap_pipeline + event_id)

# loading the MOC file
MOC_file = MOC.from_file( skymap_pipeline + event_id + '_MOC_' + str(prob_percentage) )

# square degrees in a whole sphere
from math import pi
square_degrees_sphere = (360.0**2)/pi

# printing area
area_sq2 = round( ( MOC_file.sky_fraction * square_degrees_sphere ), 1 )
print ( str( int( prob_percentage*100 ) )+'%' + ' area = ', area_sq2, 'sq. deg' )

# loading DSS colored for sky background
get_hips( "P/DSS2/color" )



3. Query Catalogs from MOCs
catalog = 'VII/267/gwgc' # selecting catalog
catalog_renamed = catalog.replace('/', '_')

# selecting MOC coverage
from mocpy import MOC
moc = MOC.from_file( 'bayestar18951_MOC_0.9' )

# sending to Aladin plane
send_file( 'bayestar18951_MOC_0.9')
rename ( 'bayestar18951_MOC_0.9' )

# querying from MOC ignoring astropy.io.votable.exceptions
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    table = moc.query_vizier_table( catalog, max_rows = 100000 ) 

# file output: votable format
table.write( catalog_renamed + 'MOC_query', format = 'votable', overwrite = True )

# sending to the Aladin plane 
send_file( catalog_renamed + 'MOC_query' )

# loading DSS colored for sky background
get_hips( "P/DSS2/color" )

3.B Ranked list of galaxies in 3D sky map
# downloading 3D HEALPix sky map
from astropy.utils.data import download_file
url = ('http://asd.gsfc.nasa.gov/Leo.Singer/'+'going-the-distance/2015/compare/18951/'+'bayestar.fits.gz')
filename = download_file(url, cache=True)

# reading HEALPix layers
import healpy as hp
prob, distmu, distsigma, distnorm = hp.read_map(filename, 
                                                field=[0, 1, 2, 3], verbose=False)

# HEALPix resolution 
npix = len(prob)
nside = hp.npix2nside(npix)

pixarea = hp.nside2pixarea(nside)

# Ranking list of galaxies from a MOC region
from mocpy import MOC
moc = MOC.from_file( 'bayestar18951_MOC_0.9' )

# sending to Aladin plane
send_file( 'bayestar18951_MOC_0.9')
rename ( 'bayestar18951_MOC_0.9' )

catalog = 'J/ApJS/199/26/table3'  # 2MASS Redshift Survey 

# querying from MOC ignoring astropy.io.votable.exceptions
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    cat = moc.query_vizier_table( catalog, max_rows = 100000 )

import numpy as np
from scipy.special import gammaincinv

completeness = 0.5
alpha = -1.0
MK_star = -23.55

MK_max = MK_star + 2.5*np.log10(gammaincinv(alpha + 2, completeness))

# selecting only galaxies with positive redshifts and absolute
# magnitudes greater than M(max)
from astropy.cosmology import WMAP9 as cosmo

import astropy.units as u
import astropy.constants as c

z = (u.Quantity(cat['cz']) / c.c).to(u.dimensionless_unscaled)

MK = cat['Ktmag'] - cosmo.distmod(z)
keep = (z > 0) & (MK < MK_max)

cat = cat[keep]
z = z[keep]

# luminosity distance and HEALPix index of each galaxy
r = cosmo.luminosity_distance(z).to('Mpc').value

theta = 0.5*np.pi - cat['_DEJ2000'].to('rad').value
phi = cat['_RAJ2000'].to('rad').value
ipix = hp.ang2pix(nside, theta, phi)

# probability density per unit volume at the position of each galaxy
from scipy.stats import norm
dp_dV = prob[ipix]*distnorm[ipix]*norm(distmu[ipix], distsigma[ipix]).pdf(r) / pixarea

#sorting the galaxies by descending probability density
galaxies_in_moc = cat[np.flipud(np.argsort(dp_dV))][:]

# adding probability galaxy position to the catalog
from astropy.table import Column

dp_dV_sort = np.flipud(np.argsort(dp_dV))[:]
dp_dV_value = dp_dV[dp_dV_sort]

# rounding
dp_dV_value_round = []
dp_dV_value_round = ['{:.3e}'.format(i) for i in dp_dV_value]

probability_galaxy_position = Column(dp_dV_value_round, name = 'dp_dV')

galaxies_in_moc.add_column(probability_galaxy_position, index=0)
print (galaxies_in_moc['_RAJ2000', '_DEJ2000', 'Ktmag','dp_dV'])

# sending to Aladin plane the weighted catalog 
galaxies_in_moc.write( 'ranked_list_galaxies', format = 'votable', overwrite = True )
send_file( 'ranked_list_galaxies' )

# loading DSS colored for sky background
get_hips( "P/DSS2/color" )

#3.C Queries running simultaneously

# selecting catalogs
catalogs = ['VII/267/gwgc','J/ApJ/675/1459/table1','VII/110A',
            'J/AJ/137/2981','J/A+A/534/A109','J/ApJS/199/26/table3'] 

# selecting MOC coverage
from mocpy import MOC
moc = MOC.from_file( 'bayestar18951_MOC_0.9') 

# sending to Aladin plane
send_file( 'bayestar18951_MOC_0.9')
rename ( 'bayestar18951_MOC_0.9' )

# querying from MOC ignoring astropy.io.votable.exceptions
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    for catalog in catalogs:
        print ()
        table = moc.query_vizier_table( catalog, max_rows = 100000 ) 
        print (table)
        catalog_renamed = catalog.replace('/', '_')
        table.write( catalog_renamed + 'MOC_query', format = 'votable', overwrite = True )
        send_file( catalog_renamed + 'MOC_query' )
        
# loading DSS colored for sky background
get_hips( "P/DSS2/color" )

## selecting MOC coverage
from mocpy import MOC
moc_1 = MOC.from_file( 'bayestar18951_MOC_0.9' ) 

# sending to Aladin plane
send_file( 'bayestar18951_MOC_0.9')
rename ( 'bayestar18951_MOC_0.9' )

# loading the MOC coverage map of SDSS Photometric Catalog (9)
from astropy.utils.data import download_file
url_id = 'http://alasky.u-strasbg.fr/footprints/tables/vizier/V_139_sdss9/MOC'
sdss9 = download_file( url_id, cache = True, timeout = 300 )
send_file( sdss9 )
rename ( 'sdss9_MOC' )

#load sdss9 MOC coverage
moc_2 = MOC.from_file( sdss9 ) 

# Intersection operation and writing file
inter = moc_1.intersection( moc_2 )
inter.write( 'inter', format = 'fits')

#sending to Aladin plane
send_file( 'inter' )
rename ( 'inter' )

# loading DSS colored for sky background
get_hips( "P/DSS2/color" )

import ipywidgets as widgets
from ipywidgets import interact, fixed

# selecting an event id (2015);
# http://www.ligo.org/scientists/first2years/
event_id = '18951'

# bayestar sky map
skymap_pipeline = 'bayestar'

# setting enclosed probability percentage 
prob_percentage = 0.9

# loading the simulated CBC event id (2015)
from astropy.utils.data import download_file

url_id = 'http://www.ligo.org/scientists/first2years/2015/compare/'+event_id+'/'+skymap_pipeline+'.fits.gz'
pipeline_event = download_file( url_id, cache = True, timeout = 300 )

# sending to Aladin plane
send_file ( pipeline_event )
rename ( skymap_pipeline + event_id )

# plotting contours from a specific url
from mocpy import MOC
get_json_link( 'https://losc.ligo.org/s/skymapViewer/json/skymaps/F2Y/'+event_id+'.json' )

#slider MOC production
interact( MOC_confidence_region, infile = pipeline_event, percentage = (0.1, 0.9, 0.1), 
         short_name = fixed( event_id ) )

# loading DSS colored for sky background
get_hips( "P/DSS2/color" )
