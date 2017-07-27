

// Fortran functions
extern "C" {
extern void initialize_(int* sz,
                        int* gc,
                        REAL_TYPE* q,
                        REAL_TYPE* w);

extern void bc_(int* sz,
                int* gc,
                REAL_TYPE* q,
                REAL_TYPE* dh,
                REAL_TYPE* org,
                int* tbl);

extern void euler_explicit_(int* sz,
                            int* gc,
                            REAL_TYPE* q,
                            REAL_TYPE* w,
                            REAL_TYPE* dh,
                            REAL_TYPE* dt,
                            REAL_TYPE* alpha,
                            REAL_TYPE* res);

extern void write_sph_(int* sz,
                       int* gc,
                       int* step,
                       REAL_TYPE* time,
                       REAL_TYPE* dh,
                       REAL_TYPE* org,
                       char* fname,
                       REAL_TYPE* q);
};
