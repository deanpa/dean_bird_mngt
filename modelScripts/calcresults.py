
import pickle

class KiwiResults(object):
    """
    Dummy class to take the parameters for the rios
    parallel processing. Needed in a separate module
    so the pickling still works.

    Contains the info for a single iteration
    """
    ## Use these results to produce mean proportion of kiwi_K within the management
    ## areas and the total area
    ## Use this to produce time series graphs for each area.
#    kiwiPropKMap_2D = None
    kiwiDensity_2D = None
    stoatDensity_2D = None
    rodentDensity_2D = None
    kiwiDensity_2D_mth = None
    stoatDensity_2D_mth = None
    rodentDensity_2D_mth = None    

    controlCount = None

    ## 3-D arrays to visual dyanamics over the years of simulation
    popAllYears_3D = None

    # parameters used for model
    params = None

    def pickleSelf(self, fname):
        fileobj = open(fname, 'wb')
        pickle.dump(self, fileobj, protocol=4) # so we get large file support
        fileobj.close()

    @staticmethod
    def unpickleFromFile(fname):
        fileobj = open(fname, 'rb')
        data = pickle.load(fileobj)
        fileobj.close()
        return data

