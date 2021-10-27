// // struct and variable definitions

// //From grids.f90

// struct _type_lonlat{
//     int * mask[:][:]; //originally logical
//     float * lonc[:][:], latc[:][:], lont[:][:], latt[:][:], vect[:][:], tic[:][:], tac[:][:], t3d[:][:];
// };

// struct _type_data{
//     int * m[:][:]; //originally logical
//     float * u[:][:], v[:][:];
// };

// struct _type_location{
//     float * ag[:][:], bg[:][:];
//     float * pg[:][:];
//     int * ig[:][:], jg[:][:];
//     int * mg[:][:]; //logical
// };

// struct _type_weights{
//     int nx, ny;
//     int numwgt;
//     float * wgt[:][:];
//     int * iwgt[:][:], jwgt[:][:];
//     int * nwgt[:][:];
//     int * mwgt[:][:]; //logical
// };


// //From kdtree2.f90

extern void cstintrp_from_call_c(double xsrc, double ysrc, double fsrc, double xdst, double ydst, double fdst);