cdef extern from "cstintrp_py.h":
    cdef void cstintrp_from_call_c(double xsrc, double ysrc, double fsrc, double xdst, double ydst, double fdst)

def cstintrp_from_call_py(double xsrc, double ysrc, double fsrc, double xdst, double ydst):
    cdef double fdst[:][:]
    
    cstintrp_from_call_c(double xsrc, double ysrc, double fsrc, double xdst, double ydst, double fdst)
