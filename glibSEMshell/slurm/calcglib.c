#pragma GCC diagnostic ignored "-Wunused-result"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <omp.h>
#include <sys/mman.h>

#define HEADER_LEN 512
#define MAX_NETWORK_LEN 8
#define MAX_STATION_LEN 32
#define STR_BUFFER_LEN 2048

#define NGLL 5
#define NDIM 3

typedef long int int64;
typedef unsigned int uint32;

#define ivec_get(ivec, idx) ivec.v[idx]

#define ivec_set(ivec, idx, v) ivec_get(ivec, idx) = v

/*
 * Int32Vec{
 *     int64 n;
 *     int *v;
 * }
 */
typedef struct int32_vector{
    int64 n;
    int* v;
} Int32Vec;

Int32Vec ivec32_alloc(int64 n){
    Int32Vec _v;
    _v.n = n;
    _v.v = (int*)malloc(n*sizeof(int));
    return _v;
}

Int32Vec ivec32_deepcopy(Int32Vec ivec){
    Int32Vec _v;
    _v.n = ivec.n;
    _v.v = (int*)malloc(ivec.n*sizeof(int));
    memcpy(_v.v, ivec.v, ivec.n*sizeof(int));
    return _v;
}

void ivec32_init(Int32Vec* ivec){
    ivec->n = 0;
    ivec->v = NULL;
    return;
}

void ivec32_free(Int32Vec* ivec){
    if(ivec->v!=NULL) free(ivec->v);
    ivec->n = 0;
    ivec->v = NULL;
    return;
}

int ivec32_prod(Int32Vec ivec){
    int v, i;
    v = 1;
    for(i = 0; i<ivec.n; i++) v *= ivec.v[i];
    return v;
}

/*
 * Int64Vec{
 *     int64 n;
 *     int64 *v;
 * }
 */
typedef struct int64_vector{
    int64 n;
    int64* v;
} Int64Vec;

Int64Vec ivec64_alloc(int64 n){
    Int64Vec _v;
    _v.n = n;
    _v.v = (int64*)malloc(n*sizeof(int64));
    return _v;
}

Int64Vec ivec64_deepcopy(Int64Vec ivec){
    Int64Vec _v;
    _v.n = ivec.n;
    _v.v = (int64*)malloc(ivec.n*sizeof(int64));
    memcpy(_v.v, ivec.v, ivec.n*sizeof(int64));
    return _v;
}

void ivec64_init(Int64Vec* ivec){
    ivec->n = 0;
    ivec->v = NULL;
    return;
}

void ivec64_free(Int64Vec* ivec){
    if(ivec->v!=NULL) free(ivec->v);
    ivec->n = 0;
    ivec->v = NULL;
    return;
}

int64 ivec64_prod(Int64Vec ivec){
    int64 v, i;
    v = 1;
    for(i = 0; i<ivec.n; i++) v *= ivec.v[i];
    return v;
}

int64 idx64_linear(Int64Vec cart, Int64Vec base){
    int64 i, idx;
    idx = cart.v[cart.n-1];
    for(i = cart.n-2; i>=0; i--){
        idx *= base.v[i];
        idx += cart.v[i];
    }
    return idx;
}

void idx64_cart(Int64Vec cart, Int64Vec base, Int64Vec page, int64 lin){
    int64 i, res;

    ivec_set(page, 0, 1);
    for(i = 1; i<base.n; i++)
        ivec_set(page, i, ivec_get(page, i-1)*ivec_get(base, i-1));
    res = lin;
    for(i = cart.n-1; i>=0; i -= 1){
        ivec_set(cart, i, res/ivec_get(page, i));
        res %= ivec_get(page, i);
    }
    return;
}

// = = = = = = = = = = = = = = = = = = = = = = = =

#define arr_get_2d(arr, i0, i1) ((arr).v[(i0)+(i1)*((arr).s.v[0])])

#define arr_get_3d(arr, i0, i1, i2) ((arr).v[(i0)+((i1)+(i2)*((arr).s.v[1]))*((arr).s.v[0])])

#define arr_get_4d(arr, i0, i1, i2, i3) ((arr).v[ \
    (i0) + ( \
        (i1) + ( \
            (i2) + ( \
                i3 \
            ) * ((arr).s.v[2]) \
        ) * ((arr).s.v[1]) \
    ) * ((arr).s.v[0]) \
])

#define arr_get_5d(arr, i0, i1, i2, i3, i4) ((arr).v[(i0)+( \
        (i1) + ( \
            (i2) + ( \
                (i3) + (i4) * ((arr).s.v[3]) \
            ) * ((arr).s.v[2]) \
        ) * ((arr).s.v[1]) \
    ) * ((arr).s.v[0])])

#define arr_get_6d(arr, i0, i1, i2, i3, i4, i5) ((arr).v[ \
    (i0) + ( \
        (i1) + ( \
            (i2) + ( \
                (i3) + ( \
                    (i4) + (i5) * ((arr).s.v[4]) \
                ) * ((arr).s.v[3]) \
            ) * ((arr).s.v[2]) \
        ) * ((arr).s.v[1]) \
    ) * ((arr).s.v[0]) \
])

/*
 * Int32Array{
 *     Int64Vec s;
 *     int *v;
 * }
 */
typedef struct i32_array{
    Int64Vec s;
    int* v;
} Int32Array;

Int32Array iarr32_alloc(Int64Vec s){
    Int32Array ia;
    ia.s = ivec64_deepcopy(s);
    ia.v = (int*)malloc(ivec64_prod(s)*sizeof(int));
    return ia;
}

void iarr32_init(Int32Array* iarr){
    ivec64_init(&(iarr->s));
    iarr->v = NULL;
    return;
}

void iarr32_read(Int32Array iarr, FILE* fp){
    fread(iarr.v, sizeof(int), ivec64_prod(iarr.s), fp);
    return;
}

void iarr32_free(Int32Array* iarr){
    ivec64_free(&(iarr->s));
    if(iarr->v!=NULL) free(iarr->v);
    iarr->v = NULL;
    return;
}

/*
 * Int64Array{
 *     Int64Vec s;
 *     int64 *v;
 * }
 */
typedef struct i64_array{
    Int64Vec s;
    int64* v;
} Int64Array;

Int64Array iarr64_alloc(Int64Vec s){
    Int64Array ia;
    ia.s = ivec64_deepcopy(s);
    ia.v = (int64*)malloc(ivec64_prod(s)*sizeof(int64));
    return ia;
}

void iarr64_init(Int64Array* iarr){
    ivec64_init(&(iarr->s));
    iarr->v = NULL;
    return;
}

void iarr64_read(Int64Array iarr, FILE* fp){
    fread(iarr.v, sizeof(int64), ivec64_prod(iarr.s), fp);
    return;
}

void iarr64_free(Int64Array* iarr){
    ivec64_free(&(iarr->s));
    if(iarr->v!=NULL) free(iarr->v);
    iarr->v = NULL;
    return;
}

/*
 * FloatArray{
 *     Int64Vec s;
 *     float *v;
 * }
 */
typedef struct f_array{
    Int64Vec s;
    float* v;
} FloatArray;

FloatArray farr_alloc(Int64Vec s){
    FloatArray fa;
    fa.s = ivec64_deepcopy(s);
    fa.v = (float*)malloc(ivec64_prod(s)*sizeof(float));
    return fa;
}

void farr_init(FloatArray* farr){
    ivec64_init(&(farr->s));
    farr->v = NULL;
    return;
}

void farr_read(FloatArray farr, FILE* fp){
    fread(farr.v, sizeof(float), ivec64_prod(farr.s), fp);
    return;
}

void farr_free(FloatArray* farr){
    ivec64_free(&(farr->s));
    if(farr->v!=NULL) free(farr->v);
    farr->v = NULL;
    return;
}

/*
 * DoubleArray{
 *     Int64Vec s;
 *     double *v;
 * }
 */
typedef struct d_array{
    Int64Vec s;
    double* v;
} DoubleArray;

DoubleArray darr_alloc(Int64Vec s){
    DoubleArray da;
    da.s = ivec64_deepcopy(s);
    da.v = (double*)malloc(ivec64_prod(s)*sizeof(double));
    return da;
}

void darr_init(DoubleArray* darr){
    ivec64_init(&(darr->s));
    darr->v = NULL;
    return;
}

void darr_read(DoubleArray darr, FILE* fp){
    fread(darr.v, sizeof(double), ivec64_prod(darr.s), fp);
    return;
}

void darr_free(DoubleArray* darr){
    ivec64_free(&(darr->s));
    if(darr->v!=NULL) free(darr->v);
    darr->v = NULL;
    return;
}

/*
 *
 * typedef struct glib_info{
 *     int id, nx, ny, nz, nt;
 *     float* x, * y, * z, * t, rt;
 * } Glib_info;
*/
typedef struct glib_info{
    int id, nx, ny, nz, nt;
    float* x, * y, * z, * t, rt;
} Glib_info;

void glib_info_init(Glib_info* g){
    g->id = -1;
    g->nx = 0;
    g->ny = 0;
    g->nz = 0;
    g->nt = 0;
    g->rt = 0.0;
    g->x = NULL;
    g->y = NULL;
    g->z = NULL;
    g->t = NULL;
    return;
}

Glib_info glib_info_read(FILE* fp){
    Glib_info g;

    glib_info_init(&g);
    fread(&(g.rt), sizeof(float), 1, fp);
    fread(&(g.nx), sizeof(int), 1, fp);
    fread(&(g.ny), sizeof(int), 1, fp);
    fread(&(g.nz), sizeof(int), 1, fp);
    fread(&(g.nt), sizeof(int), 1, fp);
    g.x = (float*)malloc(g.nx*sizeof(float));
    fread(g.x, sizeof(float), g.nx, fp);
    g.y = (float*)malloc(g.ny*sizeof(float));
    fread(g.y, sizeof(float), g.ny, fp);
    g.z = (float*)malloc(g.nz*sizeof(float));
    fread(g.z, sizeof(float), g.nz, fp);
    g.t = (float*)malloc(g.nt*sizeof(float));
    fread(g.t, sizeof(float), g.nt, fp);
    return g;
}

void glib_info_free(Glib_info* g){
    if(g->x!=NULL) free(g->x);
    if(g->y!=NULL) free(g->y);
    if(g->z!=NULL) free(g->z);
    if(g->t!=NULL) free(g->t);
    return;
}

Int32Array glib_info_read_station_locate_table(FILE* fp){
    int n_station;
    Int64Vec iv;
    Int32Array tab;

    ivec64_init(&iv);
    iarr32_init(&tab);

    fread(&n_station, sizeof(int), 1, fp);
    iv = ivec64_alloc(2);
    ivec_set(iv, 0, n_station);
    ivec_set(iv, 1, 5);
    tab = iarr32_alloc(iv);
    ivec64_free(&iv);

    iarr32_read(tab, fp);
    return tab;
}

/*
 * Header_global{
 *     int nglob;
 *     int nspec;
 *     int nt;
 *     int nrec;
 *     Int32Array spec2glob;
 *     int *rec2spec;
 *     FloatArray l;
 *     char** network;
 *     char** station;
 *     FloatArray hp;
 *     DoubleArray nu;
 * }
 */
typedef struct global_var{
    int nglob, nspec, nt, nrec;
    Int32Array spec2glob;
    int* rec2spec;
    FloatArray l;
    char** network;
    char** station;
    FloatArray hp;
    DoubleArray nu;
} Header_global;

void hdr_global_init(Header_global* hdr){
    hdr->nglob = 0;
    hdr->nspec = 0;
    hdr->nt = 0;
    hdr->nrec = 0;
    iarr32_init(&(hdr->spec2glob));
    hdr->rec2spec = NULL;
    farr_init(&(hdr->l));
    hdr->network = NULL;
    hdr->station = NULL;
    farr_init(&(hdr->hp));
    darr_init(&(hdr->nu));
    return;
}

void copy_l_part(FloatArray l, FloatArray buff, int64 i, int64 j){
    int64 i1, i2, i3, i4;
    for(i4 = 0; i4<ivec_get(l.s, 5); i4++)
        for(i3 = 0; i3<ivec_get(l.s, 4); i3++)
            for(i2 = 0; i2<ivec_get(l.s, 3); i2++)
                for(i1 = 0; i1<ivec_get(l.s, 2); i1++)
                    arr_get_6d(l, i, j, i1, i2, i3, i4) = arr_get_4d(buff, i1, i2, i3, i4);
    return;
}

void read_l(FloatArray l, FloatArray buff, FILE* fp){
    int64 i, j;
    for(j = 0; j<NDIM; j += 1)
        for(i = 0; i<NDIM; i += 1){
            farr_read(buff, fp);
            copy_l_part(l, buff, i, j);
        }
    return;
}

void hdr_global_read(Header_global* hdr, FILE* fp){
    int i_buff[8], i, j, k;
    float f_buff;
    Int64Vec iv;
    FloatArray fa;

    ivec64_init(&iv);
    farr_init(&fa);
    // 8 int
    fread(i_buff, sizeof(int), 8, fp);
    if(i_buff[0]!=NDIM){
        printf("(hdr_global_read) Dimension error");
        exit(0);
    }
    hdr->nglob = i_buff[1];
    hdr->nspec = i_buff[2];
    hdr->nt = i_buff[3];
    if((i_buff[4]!=NGLL)||(i_buff[5]!=NGLL)||(i_buff[6]!=NGLL)){
        printf("NGLL error");
        exit(0);
    }
    hdr->nrec = i_buff[7];
    // spec2glob
    iv = ivec64_alloc(4);
    ivec_set(iv, 0, NGLL);
    ivec_set(iv, 1, NGLL);
    ivec_set(iv, 2, NGLL);
    ivec_set(iv, 3, hdr->nspec);
    hdr->spec2glob = iarr32_alloc(iv);
    ivec64_free(&iv);
    iarr32_read(hdr->spec2glob, fp);
    // rec2spec
    hdr->rec2spec = (int*)malloc((hdr->nrec)*sizeof(int));
    fread(hdr->rec2spec, sizeof(int), hdr->nrec, fp);
    // l
    iv = ivec64_alloc(6);
    ivec_set(iv, 0, NDIM);
    ivec_set(iv, 1, NDIM);
    ivec_set(iv, 2, NGLL);
    ivec_set(iv, 3, NGLL);
    ivec_set(iv, 4, NGLL);
    ivec_set(iv, 5, hdr->nspec);
    hdr->l = farr_alloc(iv);
    ivec64_free(&iv);

    iv = ivec64_alloc(4);
    ivec_set(iv, 0, NGLL);
    ivec_set(iv, 1, NGLL);
    ivec_set(iv, 2, NGLL);
    ivec_set(iv, 3, hdr->nspec);
    fa = farr_alloc(iv);
    ivec64_free(&iv);

    read_l(hdr->l, fa, fp);
    // network
    hdr->network = (char**)malloc((hdr->nrec)*sizeof(char*));
    for(i = 0; i<hdr->nrec; i += 1){
        (hdr->network)[i] = (char*)malloc((MAX_NETWORK_LEN+1)*sizeof(char));
        fread((hdr->network)[i], sizeof(char), MAX_NETWORK_LEN, fp);
        (hdr->network)[i][MAX_NETWORK_LEN] = '\0';
    }
    // station
    hdr->station = (char**)malloc((hdr->nrec)*sizeof(char*));
    for(i = 0; i<hdr->nrec; i += 1){
        (hdr->station)[i] = (char*)malloc((MAX_STATION_LEN+1)*sizeof(char));
        fread((hdr->station)[i], sizeof(char), MAX_STATION_LEN, fp);
        (hdr->station)[i][MAX_STATION_LEN] = '\0';
    }
    // hp
    iv = ivec64_alloc(3);
    ivec_set(iv, 0, NDIM);
    ivec_set(iv, 1, NGLL);
    ivec_set(iv, 2, NGLL);
    hdr->hp = farr_alloc(iv);
    for(i = 0; i<NDIM; i++)
        for(j = 0; j<NGLL; j++)
            for(k = 0; k<NGLL; k++)
                fread(&(arr_get_3d(hdr->hp, i, j, k)), sizeof(float), 1, fp);
    ivec64_free(&iv);
    // nu

    iv = ivec64_alloc(3);
    ivec_set(iv, 0, NDIM);
    ivec_set(iv, 1, NDIM);
    ivec_set(iv, 2, hdr->nrec);
    hdr->nu = darr_alloc(iv);
    ivec64_free(&iv);
    darr_read(hdr->nu, fp);
    return;
}

void hdr_global_free(Header_global* hdr){
    int i;
    iarr32_free(&(hdr->spec2glob));
    farr_free(&(hdr->l));
    farr_free(&(hdr->hp));
    darr_free(&(hdr->nu));
    if(hdr->rec2spec!=NULL) free(hdr->rec2spec);
    hdr->rec2spec = NULL;
    for(i = 0; i<hdr->nrec; i++){
        if((hdr->network)[i]!=NULL) free((hdr->network)[i]);
        if((hdr->station)[i]!=NULL) free((hdr->station)[i]);
    }
    if(hdr->network!=NULL) free(hdr->network);
    hdr->network = NULL;
    if(hdr->station!=NULL) free(hdr->station);
    hdr->station = NULL;
    hdr->nglob = 0;
    hdr->nspec = 0;
    hdr->nt = 0;
    hdr->nrec = 0;
    return;
}

const char glibio_digits[62] = { '0', '1', '2', '3', '4', '5', '6', '7',
                                '8', '9', 'a', 'b', 'c', 'd', 'e', 'f',
                                'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
                                'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                                'w', 'x', 'y', 'z', 'A', 'B', 'C', 'D',
                                'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L',
                                'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',
                                'U', 'V', 'W', 'X', 'Y', 'Z' };
const int64 glibio_K = 62;

Int64Vec hdr_decode_station_id(Header_global ghdr){
    int64 i, j, k, l, v;
    Int64Vec iv;
    char strbuff[MAX_NETWORK_LEN+MAX_STATION_LEN+1] = { '\0' };
    ivec64_init(&iv);
    iv = ivec64_alloc((int64)ghdr.nrec);
#pragma omp parallel for default(none) private(i, j, k, l, strbuff, v) shared(ghdr, iv)
    for(i = 0; i<ghdr.nrec; i++){
        for(j = 0; j<MAX_NETWORK_LEN; j++) strbuff[j] = ghdr.network[i][j];
        for(j = 0; j<MAX_STATION_LEN; j++) strbuff[j+MAX_NETWORK_LEN] = ghdr.station[i][j];
        for(j = 0; j<MAX_NETWORK_LEN+MAX_STATION_LEN; j++) if(strbuff[j]!='0') break;
        v = 0;
        for(k = 0; k<MAX_NETWORK_LEN+MAX_STATION_LEN; k++){
            v *= glibio_K;
            for(l = 0; l<glibio_K; l++) if(glibio_digits[l]==strbuff[k]) break;
            v += l;
        }
        ivec_set(iv, i, v);
    }
    return iv;
}

/*
 * Header_local{
 *     int lrec;
 *     int lglob;
 *     int* lglob2glob;
 *     int* lrec2rec;
 *     FloatArray xi_r;
 *     FloatArray eta_r;
 *     FloatArray gamma_r;
 * }
 */
typedef struct header_local{
    int lrec, lglob;
    int* lglob2glob;
    int* lrec2rec;
    FloatArray xi_r, eta_r, gamma_r;
} Header_local;

void hdr_local_init(Header_local* hdr){
    hdr->lrec = 0;
    hdr->lglob = 0;
    hdr->lglob2glob = NULL;
    hdr->lrec2rec = NULL;
    farr_init(&(hdr->xi_r));
    farr_init(&(hdr->eta_r));
    farr_init(&(hdr->gamma_r));
    return;
}

void hdr_local_read(Header_local* hdr, FILE* fp){
    int i_buff[5];
    Int64Vec iv;

    ivec64_init(&iv);
    fread(i_buff, sizeof(int), 5, fp);
    if(i_buff[0]!=NDIM){
        printf("Dimension error");
        exit(0);
    }
    if(i_buff[1]!=NGLL){
        printf("NGLL error");
        exit(0);
    }
    hdr->lrec = i_buff[3];
    hdr->lglob = i_buff[4];
    hdr->lglob2glob = (int*)malloc(hdr->lglob*sizeof(int));
    fread(hdr->lglob2glob, sizeof(int), hdr->lglob, fp);
    hdr->lrec2rec = (int*)malloc(hdr->lrec*sizeof(int));
    fread(hdr->lrec2rec, sizeof(int), hdr->lrec, fp);

    iv = ivec64_alloc(2);
    ivec_set(iv, 0, NGLL);
    ivec_set(iv, 1, hdr->lrec);
    hdr->xi_r = farr_alloc(iv);
    farr_read(hdr->xi_r, fp);
    hdr->eta_r = farr_alloc(iv);
    farr_read(hdr->eta_r, fp);
    hdr->gamma_r = farr_alloc(iv);
    farr_read(hdr->gamma_r, fp);
    ivec64_free(&iv);
    return;
}

void hdr_local_free(Header_local* hdr){
    farr_free(&(hdr->xi_r));
    farr_free(&(hdr->eta_r));
    farr_free(&(hdr->gamma_r));
    if(hdr->lglob2glob!=NULL) free(hdr->lglob2glob);
    hdr->lglob2glob = NULL;
    if(hdr->lrec2rec!=NULL) free(hdr->lrec2rec);
    hdr->lrec2rec = NULL;
    hdr->lrec = 0;
    hdr->lglob = 0;
    return;
}

FloatArray local_field_read(FILE* fp){
    int i_buff[3];
    Int64Vec iv;
    FloatArray _fa;

    ivec64_init(&iv);

    fread(i_buff, sizeof(int), 3, fp);
    if(i_buff[0]!=NDIM){
        printf("Dimension error");
        exit(0);
    }
    iv = ivec64_alloc(3);
    ivec_set(iv, 0, i_buff[0]);
    ivec_set(iv, 1, i_buff[2]);
    ivec_set(iv, 2, i_buff[1]);
    _fa = farr_alloc(iv);
    ivec64_free(&iv);
    farr_read(_fa, fp);
    return _fa;
}


typedef struct interp_buffer{
    FloatArray b0, b1, b2, b3, duidxj;
    int nglob, lglob, nspec, lspec, nrec, lrec, nt, nthreads;
    int* glob2lglob, * lspec2spec, * spec2lspec;
} Interp_buffer;

void interp_init_buffer(Interp_buffer* buffer){
    buffer->nglob = 0;
    buffer->lglob = 0;
    buffer->nspec = 0;
    buffer->lspec = 0;
    buffer->nrec = 0;
    buffer->lrec = 0;
    buffer->nt = 0;
    buffer->nthreads = 0;
    buffer->glob2lglob = NULL;
    buffer->spec2lspec = NULL;
    buffer->lspec2spec = NULL;
    farr_init(&(buffer->b0));
    farr_init(&(buffer->b1));
    farr_init(&(buffer->b2));
    farr_init(&(buffer->b3));
    farr_init(&(buffer->duidxj));
    return;
}

void interp_free_buffer(Interp_buffer* buffer){
    int64 i;

    farr_free(&(buffer->b0));
    farr_free(&(buffer->b1));
    farr_free(&(buffer->b2));
    farr_free(&(buffer->b3));
    farr_free(&(buffer->duidxj));
    if(buffer->glob2lglob!=NULL) free(buffer->glob2lglob);
    if(buffer->spec2lspec!=NULL) free(buffer->spec2lspec);
    if(buffer->lspec2spec!=NULL) free(buffer->lspec2spec);
    buffer->nglob = 0;
    buffer->lglob = 0;
    buffer->nspec = 0;
    buffer->lspec = 0;
    buffer->nrec = 0;
    buffer->lrec = 0;
    buffer->nt = 0;
    buffer->nthreads = 0;
    return;
}

void interp_allocate_buffer(FloatArray* green, Interp_buffer* buffer, Header_global ghdr, Header_local lhdr,
    FloatArray field){
    int irl, ir, ispec, ispecl, iglob, iglobl, nglob, lglob, nspec, lspec, nrec, lrec, nt, nthread, ithread;
    Int64Vec iv;
    ivec64_init(&iv);

    nglob = ghdr.nglob;
    lglob = lhdr.lglob;
    nspec = ghdr.nspec;
    nrec = ghdr.nrec;
    lrec = lhdr.lrec;
    nt = ivec_get(field.s, 2);
    nthread = omp_get_max_threads();

    buffer->glob2lglob = (int*)malloc(nglob*sizeof(int));
    for(iglobl = 0; iglobl<lglob; iglobl++){
        iglob = lhdr.lglob2glob[iglobl];
        (buffer->glob2lglob)[iglob-1] = iglobl+1;
    }

    iv = ivec64_alloc(nspec);
    for(ispec = 0; ispec<nspec; ispec++) ivec_set(iv, ispec, 0);
    for(irl = 0; irl<lrec; irl++){
        ir = lhdr.lrec2rec[irl]-1;
        ispec = ghdr.rec2spec[ir]-1;
        ivec_set(iv, ispec, 1);
    }

    lspec = 0;
    for(ispec = 0; ispec<nspec; ispec++) lspec += ivec_get(iv, ispec);

    buffer->lspec2spec = (int*)malloc(lspec*sizeof(int));
    buffer->spec2lspec = (int*)malloc(nspec*sizeof(int));
    memset(buffer->lspec2spec, 0, lspec*sizeof(int));
    memset(buffer->spec2lspec, 0, nspec*sizeof(int));
    ispecl = 0;
    for(ispec = 0; ispec<nspec; ispec++){
        if(ivec_get(iv, ispec)>0){
            (buffer->lspec2spec)[ispecl] = ispec+1;
            (buffer->spec2lspec)[ispec] = ispecl+1;
            ispecl += 1;
        }
    }
    ivec64_free(&iv);

    buffer->nglob = nglob;
    buffer->lglob = lglob;
    buffer->nspec = nspec;
    buffer->lspec = lspec;
    buffer->nrec = nrec;
    buffer->lrec = lrec;
    buffer->nt = nt;
    buffer->nthreads = nthread;

    iv = ivec64_alloc(2);
    ivec_set(iv, 0, NDIM);
    ivec_set(iv, 1, lhdr.lrec);
    buffer->b0 = farr_alloc(iv);
    ivec64_free(&iv);

    iv = ivec64_alloc(3);
    ivec_set(iv, 0, NDIM);
    ivec_set(iv, 1, NDIM);
    ivec_set(iv, 2, lrec);
    buffer->duidxj = farr_alloc(iv);

    ivec_set(iv, 0, nt);
    ivec_set(iv, 1, 6);
    ivec_set(iv, 2, lhdr.lrec);
    *green = farr_alloc(iv);
    ivec64_free(&iv);

    iv = ivec64_alloc(4);
    ivec_set(iv, 0, NGLL);
    ivec_set(iv, 1, NGLL);
    ivec_set(iv, 2, NGLL);
    ivec_set(iv, 3, lspec);
    buffer->b1 = farr_alloc(iv);
    ivec64_free(&iv);

    iv = ivec64_alloc(5);
    ivec_set(iv, 0, NDIM);
    ivec_set(iv, 1, NGLL);
    ivec_set(iv, 2, NGLL);
    ivec_set(iv, 3, NGLL);
    ivec_set(iv, 4, lspec);
    buffer->b2 = farr_alloc(iv);
    buffer->b3 = farr_alloc(iv);
    ivec64_free(&iv);
    return;
}

void
interp_interpolate(FloatArray green, Header_global ghdr, Header_local lhdr, FloatArray field, Interp_buffer buffer){
    int i_glob, i_glob_l, i_spec, i_spec_l, ir, irl, i_time;
    int i, j, k, l, i_dim, j_dim, c_dim;
    float f_buff;
    double d_buff;
    for(i_time = 0; i_time<buffer.nt; i_time++){
        for(c_dim = 0; c_dim<NDIM; c_dim++){

#pragma omp parallel for default(none) private(i_glob, i_glob_l, i_spec, i_spec_l, i, j, k, l, i_dim, j_dim, f_buff) \
    shared(buffer, ghdr, lhdr, field, c_dim, i_time)
            for(i_spec_l = 0; i_spec_l<buffer.lspec; i_spec_l++){
                i_spec = buffer.lspec2spec[i_spec_l]-1;

                for(k = 0; k<NGLL; k++)
                    for(j = 0; j<NGLL; j++)
                        for(i = 0; i<NGLL; i++){
                            i_glob = arr_get_4d(ghdr.spec2glob, i, j, k, i_spec)-1;
                            i_glob_l = buffer.glob2lglob[i_glob]-1;
                            arr_get_4d(buffer.b1, i, j, k, i_spec_l) = arr_get_3d(field, c_dim, i_glob_l, i_time);
                        }

                for(k = 0; k<NGLL; k++)
                    for(j = 0; j<NGLL; j++)
                        for(i = 0; i<NGLL; i++)
                            for(l = 0; l<NDIM; l++){
                                arr_get_5d(buffer.b2, l, i, j, k, i_spec_l) = 0.0;
                                arr_get_5d(buffer.b3, l, i, j, k, i_spec_l) = 0.0;
                            }

                for(k = 0; k<NGLL; k++)
                    for(j = 0; j<NGLL; j++)
                        for(i = 0; i<NGLL; i++)
                            for(l = 0; l<NGLL; l++){
                                arr_get_5d(buffer.b2, 0, i, j, k, i_spec_l) += arr_get_4d(buffer.b1, l, j, k, i_spec_l)*arr_get_3d(ghdr.hp, 0, l, i);
                                arr_get_5d(buffer.b2, 1, i, j, k, i_spec_l) += arr_get_4d(buffer.b1, i, l, k, i_spec_l)*arr_get_3d(ghdr.hp, 1, l, j);
                                arr_get_5d(buffer.b2, 2, i, j, k, i_spec_l) += arr_get_4d(buffer.b1, i, j, l, i_spec_l)*arr_get_3d(ghdr.hp, 2, l, k);
                            }

                for(k = 0; k<NGLL; k++)
                    for(j = 0; j<NGLL; j++)
                        for(i = 0; i<NGLL; i++)
                            for(i_dim = 0; i_dim<NDIM; i_dim++)
                                for(j_dim = 0; j_dim<NDIM; j_dim++){
                                    arr_get_5d(buffer.b3, i_dim, i, j, k, i_spec_l) +=
                                        arr_get_5d(buffer.b2, j_dim, i, j, k, i_spec_l)*
                                        arr_get_6d(ghdr.l, i_dim, j_dim, i, j, k, i_spec);
                                }
            }

#pragma omp parallel for default(none) private(ir, irl, i_spec, i_spec_l, i_dim, j_dim, i, j, k, f_buff, d_buff) \
            shared(lhdr, ghdr, buffer, c_dim)
            for(irl = 0; irl<lhdr.lrec; irl++){
                ir = lhdr.lrec2rec[irl]-1;
                i_spec = ghdr.rec2spec[ir]-1;
                i_spec_l = buffer.spec2lspec[i_spec]-1;
                for(i_dim = 0; i_dim<NDIM; i_dim++){
                    arr_get_2d(buffer.b0, i_dim, irl) = 0.0;
                    for(k = 0; k<NGLL; k++)
                        for(j = 0; j<NGLL; j++)
                            for(i = 0; i<NGLL; i++){
                                arr_get_2d(buffer.b0, i_dim, irl) +=
                                    arr_get_5d(buffer.b3, i_dim, i, j, k, i_spec_l)*
                                    arr_get_2d(lhdr.xi_r, i, irl)*
                                    arr_get_2d(lhdr.eta_r, j, irl)*
                                    arr_get_2d(lhdr.gamma_r, k, irl);
                            }
                }

                for(i_dim = 0; i_dim<NDIM; i_dim++){
                    // arr_get_3d(buffer.duidxj, c_dim, i_dim, irl) = 0.0;
                    d_buff = 0.0;
                    for(j_dim = 0; j_dim<NDIM; j_dim++){
                        d_buff += arr_get_3d(ghdr.nu, i_dim, j_dim, ir)*arr_get_2d(buffer.b0, j_dim, irl);
                    }
                    arr_get_3d(buffer.duidxj, c_dim, i_dim, irl) = d_buff;
                }
            }
        }
#pragma omp parallel for default(none) private(irl, ir, f_buff) shared(lhdr, buffer, i_time, green)
        for(irl = 0; irl<lhdr.lrec; irl++){
            arr_get_3d(green, i_time, 0, irl) = arr_get_3d(buffer.duidxj, 1, 1, irl);
            arr_get_3d(green, i_time, 1, irl) = arr_get_3d(buffer.duidxj, 0, 0, irl);
            arr_get_3d(green, i_time, 2, irl) = arr_get_3d(buffer.duidxj, 2, 2, irl);
            arr_get_3d(green, i_time, 3, irl) = arr_get_3d(buffer.duidxj, 0, 1, irl)+arr_get_3d(buffer.duidxj, 1, 0, irl);
            arr_get_3d(green, i_time, 4, irl) = -(arr_get_3d(buffer.duidxj, 2, 1, irl)+arr_get_3d(buffer.duidxj, 1, 2, irl));
            arr_get_3d(green, i_time, 5, irl) = -(arr_get_3d(buffer.duidxj, 2, 0, irl)+arr_get_3d(buffer.duidxj, 0, 2, irl));
        }
    }

}

void write_green_to_file(char* root, FloatArray green, Int64Vec station_id_list, Int32Array station_table,
    Header_global ghdr, Header_local lhdr, int maxglib, int icmp){
    float rt, * gbuff;
    int64 shift, nlen;
    unsigned long int buffsize;
    int i, irl, ir, rid, tid, iglib, used_glib, nx, ny, nz, nt, nt0, ix, iy, iz;
    int* glib_ids, * table_ids, * unique_ids, ** id_set, * glib_counts, * current_irec;
    char strbuf[STR_BUFFER_LEN];
    FILE* fp;
    int fd;


    glib_ids = (int*)malloc(lhdr.lrec*sizeof(int));
    table_ids = (int*)malloc(lhdr.lrec*sizeof(int));
    glib_counts = (int*)malloc(maxglib*sizeof(int));
    unique_ids = (int*)malloc(maxglib*sizeof(int));

    used_glib = 0;
    for(iglib = 0; iglib<maxglib; iglib++) glib_counts[iglib] = 0;
    for(irl = 0; irl<lhdr.lrec; irl++){
        ir = lhdr.lrec2rec[irl]-1;
        rid = ivec_get(station_id_list, ir);
        table_ids[irl] = rid-1;
        glib_ids[irl] = arr_get_2d(station_table, rid-1, 1);
        for(iglib = 0; iglib<used_glib; iglib++) if(glib_ids[irl]==unique_ids[iglib]) break;
        if(iglib==used_glib){
            unique_ids[used_glib] = glib_ids[irl];
            used_glib += 1;
        }
        glib_counts[iglib] += 1;
    }

    id_set = (int**)malloc(used_glib*sizeof(int*));
    current_irec = (int*)malloc(used_glib*sizeof(int));
    for(iglib = 0; iglib<used_glib; iglib++){
        id_set[iglib]==NULL;
        id_set[iglib] = (int*)malloc(glib_counts[iglib]*sizeof(int));
        if(id_set[iglib]==NULL){
            printf("alloc id_set error\n");
            exit(0);
        }
        current_irec[iglib] = 0;
    }
    for(irl = 0; irl<lhdr.lrec; irl++)
        for(iglib = 0; iglib<used_glib; iglib++)
            if(glib_ids[irl]==unique_ids[iglib]){
                id_set[iglib][current_irec[iglib]] = irl;
                current_irec[iglib] += 1;
                break;
            }

    nt = ivec_get(green.s, 0);
    nlen = nt*6;
    for(iglib = 0; iglib<used_glib; iglib++){
        sprintf(strbuf, "%s/glib_tmp2_%d.bin", root, unique_ids[iglib]);
        fp = fopen(strbuf, "r");
        fread(&rt, sizeof(float), 1, fp);
        fread(&nx, sizeof(int), 1, fp);
        fread(&ny, sizeof(int), 1, fp);
        fread(&nz, sizeof(int), 1, fp);
        fread(&nt0, sizeof(int), 1, fp);
        fclose(fp);
        if(nt!=nt0){
            printf("nt not match, nt: %d, nt0: %d\n", nt, nt0);
            exit(0);
        }
        buffsize = (unsigned long int)nx*ny*nz*nt*18;
        sprintf(strbuf, "%s/glib_tmp3_%d.bin", root, unique_ids[iglib]);
        fd = open(strbuf, O_RDWR);
        gbuff = (float*)mmap(NULL, buffsize*sizeof(float), PROT_READ|PROT_WRITE, MAP_SHARED, fd, 0);
        for(i = 0; i<glib_counts[iglib]; i++){
            irl = id_set[iglib][i];
            iz = arr_get_2d(station_table, table_ids[irl], 2)-1;
            iy = arr_get_2d(station_table, table_ids[irl], 3)-1;
            ix = arr_get_2d(station_table, table_ids[irl], 4)-1;
            shift = (icmp+(iz+(iy+ix*(int64)ny)*nz)*3)*6*nt;
            memcpy(&(gbuff[shift]), &(arr_get_3d(green, 0, 0, irl)), nlen*sizeof(float));
        }
        munmap(gbuff, buffsize*sizeof(float));
        close(fd);
    }

    for(iglib = 0; iglib<used_glib; iglib++) free(id_set[iglib]);
    free(glib_ids);
    free(table_ids);
    free(glib_counts);
    free(unique_ids);
    free(id_set);
    free(current_irec);
    return;
}

void printhelp(void){
    printf("usage:\n");
    printf("\tcalc_glib path nproc nglib nthread\n");
    return;
}

/*
 * root nproc nglib nthread
 */
int main(int argc, char* argv[]){
    int n_thread, n_proc, n_glib, iglib, iproc, irl, ir, id, ix, iy, iz, icmp;
    int64 nlen, nshift;
    char root_path[STR_BUFFER_LEN] = "/data/testsemshell/compress",
        path_buffer[STR_BUFFER_LEN], cmp[3] = { 'N', 'E', 'D' };
    Glib_info* glib_headers;
    Header_global ghdr;
    Header_local lhdr;
    FloatArray field, green;
    Interp_buffer itp_buffer;
    Int32Array station_table;
    Int64Vec station_id_list, * sizes = NULL, coor, iv2, iv3;
    FILE* fp_in;
    double timesum[4] = { 0.0 }, t1, t2;

    if((argc>5)||(argc<4)){
        printhelp();
        exit(0);
    }
    if(argc<5)
        n_thread = 4;
    else
        n_thread = atoi(argv[4]);
    n_glib = atoi(argv[3]);
    n_proc = atoi(argv[2]);
    strcpy(root_path, argv[1]);

    printf("root dir: %s\n", root_path);
    printf("nproc: %d, nglib: %d, nthread: %d\n", n_proc, n_glib, n_thread);

    printf("init\n");
    omp_set_num_threads(n_thread);

    hdr_global_init(&ghdr);
    hdr_local_init(&lhdr);
    interp_init_buffer(&itp_buffer);
    farr_init(&field);
    farr_init(&green);
    iarr32_init(&station_table);
    ivec64_init(&station_id_list);
    ivec64_init(&coor);
    ivec64_init(&iv2);
    ivec64_init(&iv3);

    printf("read station table\n");
    sprintf(path_buffer, "%s/station_id_table.bin", root_path);
    fp_in = fopen(path_buffer, "r");
    station_table = glib_info_read_station_locate_table(fp_in);
    fclose(fp_in);

    printf("open glib file\n");

    glib_headers = (Glib_info*)malloc(n_glib*sizeof(Glib_info));
    sizes = (Int64Vec*)malloc(n_glib*sizeof(Int64Vec));
    for(iglib = 0; iglib<n_glib; iglib++){
        sprintf(path_buffer, "%s/glib_tmp2_%d.bin", root_path, iglib+1);
        fp_in = fopen(path_buffer, "r");
        glib_headers[iglib] = glib_info_read(fp_in);
        fclose(fp_in);

        nlen = glib_headers[iglib].nx*
            glib_headers[iglib].ny*
            glib_headers[iglib].nz*
            glib_headers[iglib].nt*18;
        sizes[iglib] = ivec64_alloc(6);
        ivec_set(sizes[iglib], 0, glib_headers[iglib].nt);
        ivec_set(sizes[iglib], 1, 6);
        ivec_set(sizes[iglib], 2, 3);
        ivec_set(sizes[iglib], 3, glib_headers[iglib].nz);
        ivec_set(sizes[iglib], 4, glib_headers[iglib].ny);
        ivec_set(sizes[iglib], 5, glib_headers[iglib].nx);

        field = farr_alloc(sizes[iglib]);
        sprintf(path_buffer, "%s/glib_tmp3_%d.bin", root_path, iglib+1);
        fp_in = fopen(path_buffer, "w");
        fwrite(field.v, sizeof(float), ivec64_prod(field.s), fp_in);
        fclose(fp_in);
        farr_free(&field);
    }

    printf("read global info\n");
    sprintf(path_buffer, "%s/shotT/OUTPUT_FILES/globinterp.bin", root_path);
    fp_in = fopen(path_buffer, "r");
    hdr_global_read(&ghdr, fp_in);
    fclose(fp_in);
    station_id_list = hdr_decode_station_id(ghdr);

    printf("loop on proc\n");
    iv2 = ivec64_alloc(2);
    iv3 = ivec64_alloc(3);
    coor = ivec64_alloc(6);
    for(iproc = 0; iproc<n_proc; iproc++){
        t1 = omp_get_wtime();
        sprintf(path_buffer, "%s/shotT/OUTPUT_FILES/proc%06dlocalinterp.bin", root_path, iproc);
        fp_in = fopen(path_buffer, "r");
        hdr_local_read(&lhdr, fp_in);
        fclose(fp_in);
        if(lhdr.lrec==0) continue;
        t2 = omp_get_wtime();
        timesum[0] += t2-t1;

        for(icmp = 0; icmp<3; icmp++){
            t1 = omp_get_wtime();
            printf("  cmp %c proc %d of %d\n", cmp[icmp], iproc, n_proc);
            // shot D, E, N
            sprintf(path_buffer, "%s/shot%c/OUTPUT_FILES/proc%06ddispl.bin", root_path, cmp[icmp], iproc);

            fp_in = fopen(path_buffer, "r");
            field = local_field_read(fp_in);
            fclose(fp_in);
            t2 = omp_get_wtime();
            timesum[0] += t2-t1;
            // printf("    time used: %fms\n", (t2-t1)*1000.0);

            t1 = omp_get_wtime();
            interp_allocate_buffer(&green, &itp_buffer, ghdr, lhdr, field);
            t2 = omp_get_wtime();
            timesum[1] += t2-t1;
            // printf("    time used: %fms\n", (t2-t1)*1000.0);

            t1 = omp_get_wtime();
            interp_interpolate(green, ghdr, lhdr, field, itp_buffer);
            t2 = omp_get_wtime();
            timesum[2] += t2-t1;
            // printf("    time used: %fms\n", (t2-t1)*1000.0);
            farr_free(&field);
            interp_free_buffer(&itp_buffer);

            t1 = omp_get_wtime();
            write_green_to_file(root_path, green, station_id_list, station_table, ghdr, lhdr,
                n_glib, icmp);
            t2 = omp_get_wtime();
            timesum[3] += t2-t1;
            // printf("    time used: %fms\n", (t2-t1)*1000.0);
            farr_free(&green);
        }
        hdr_local_free(&lhdr);
    }

    printf("close file and free memory\n");
    ivec64_free(&coor);
    ivec64_free(&iv2);
    ivec64_free(&iv3);
    for(iglib = 0; iglib<n_glib; iglib++){
        glib_info_free(&(glib_headers[iglib]));
        ivec64_free(&(sizes[iglib]));
    }
    free(glib_headers);
    free(sizes);
    hdr_global_free(&ghdr);
    iarr32_free(&station_table);
    ivec64_free(&station_id_list);
    printf("finish\n");

    t1 = 0;
    for(ix = 0; ix<4; ix++) t1 += timesum[ix];
    printf("time percent read:   %f%%\n", timesum[0]*100.0/t1);
    printf("time percent alloc:  %f%%\n", timesum[1]*100.0/t1);
    printf("time percent interp: %f%%\n", timesum[2]*100.0/t1);
    printf("time percent write:  %f%%\n", timesum[3]*100.0/t1);

    return 0;
}
