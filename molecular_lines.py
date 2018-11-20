from ctypes import CDLL, c_void_p, c_char_p, c_int, c_float, c_double, \
                   c_ulonglong, create_string_buffer, byref
from numpy.ctypeslib import ndpointer


#These must stay in-sync with the c library.
ml = CDLL("./libmolecular_lines.so")
fp_t = c_double
H2O = 1
CO2 = 2
O3 = 3
N2O = 4
CO = 5
CH4 = 6
O2 = 7


def get_ctype(val,t,pointer=False):
    """
    Utility function to handle arguments that are "optional".
    """
    if val == None:
        return val
    if t == "string":
        v = c_char_p(str.encode(val))
    elif t == "int":
        v = c_int(val)
    elif t == "float":
        v = c_float(val)
    elif t == "double":
        v = c_double(val)
    else:
        raise ValueError("Unsupported type.")
    if pointer:
        v = byref(v)
    return v


class MolecularLineError(BaseException):
    """
    Base exception for all libmolecular_lines errors.
    """
    pass


def catch(rc):
    """
    Check the return code from the c functions and raise the appropriate
    exception.
    """
    if rc != 0:
        s = 256
        b = create_string_buffer(s)
        ml.grt_errstr(rc,b,s)
        raise MolecularLineError(b.value.decode())


class MolecularLines(object):

    def __init__(self, num_levels, w0, wn, wres, hitran_path,
                 h2o_ctm_dir=None, o3_ctm_dir=None, wcutoff=None,
                 gpu_id=None, num_threads=None, optical_depth_method=None):
        """
        Construct the object.  The address of memory that is reserved on the
        c-side is stored in self.context.  Losing this address creates a
        memory leak.
        """
        self.ml = ml
        self.context = c_void_p(None)
        h2o_ctm_dir = get_ctype(h2o_ctm_dir,"string")
        o3_ctm_dir = get_ctype(o3_ctm_dir,"string")
        wcutoff = get_ctype(wcutoff,"double",pointer=True)
        gpu_id = get_ctype(gpu_id,"int",pointer=True)
        num_threads = get_ctype(num_threads,"int",pointer=True)
        optical_depth_method = get_ctype(optical_depth_method,"int",
                                         pointer=True)
        catch(self.ml.grt_context_init(byref(self.context),
                                       c_int(num_levels),
                                       c_double(w0),
                                       c_double(wn),
                                       c_double(wres),
                                       c_char_p(str.encode(hitran_path)),
                                       h2o_ctm_dir,
                                       o3_ctm_dir,
                                       wcutoff,
                                       gpu_id,
                                       num_threads,
                                       optical_depth_method))

    def __del__(self):
        """
        Release memory that was malloc'd on the c-side.
        """
        self.ml.grt_context_free(byref(self.context))

    def add_molecule(self, molecule_id, min_line_center=None,
                     max_line_center=None):
        """
        Wrapper for grt_add_molecule.
        """
        min_line_center = get_ctype(min_line_center,"double",pointer=True)
        max_line_center = get_ctype(max_line_center,"double",pointer=True)
        catch(self.ml.grt_add_molecule(self.context,
                                       c_int(molecule_id),
                                       min_line_center,
                                       max_line_center))

    def set_ppmv(self, molecule_id, ppmv):
        """
        Wrapper for grt_set_molecule_ppmv.
        """
        self.ml.grt_set_molecule_ppmv.argtypes = [c_void_p,
                                                  c_int,
                                                  ndpointer(fp_t, flags="C_CONTIGUOUS")]
        catch(self.ml.grt_set_molecule_ppmv(self.context,
                                            c_int(molecule_id),
                                            ppmv))

    def optical_depths(self, pressure, temperature, optical_depth):
        """
        Wrapper for grt_calculate_optical_depth.
        """
        self.ml.grt_calculate_optical_depth.argtypes = [c_void_p,
                                                        ndpointer(fp_t, flags="C_CONTIGUOUS"),
                                                        ndpointer(fp_t, flags="C_CONTIGUOUS"),
                                                        ndpointer(fp_t, flags="C_CONTIGUOUS")]
        catch(self.ml.grt_calculate_optical_depth(self.context,
                                                  pressure,
                                                  temperature,
                                                  optical_depth))

    def num_molecules(self):
        """
        Wrapper for grt_get_num_molecules.
        """
        n = c_int()
        catch(self.ml.grt_get_num_molecules(self.context,
                                            byref(n)))
        return n.value

    def spectral_grid_size(self):
        """
        Wrapper for grt_get_spectral_grid_size.
        """
        n = c_ulonglong()
        catch(self.ml.grt_get_spectral_grid_size(self.context,
                                                 byref(n)))
        return n.value


def set_verbosity(level):
    """
    Wrapper for grt_set_verbosity.
    """
    n = c_int(level)
    ml.grt_set_verbosity(n)


def verbosity():
    """
    Wrapper for grt_get_verbosity.
    """
    return ml.grt_get_verbosity()
