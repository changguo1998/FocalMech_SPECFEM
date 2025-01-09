#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "omp.h"

typedef __uint8_t uint8;
typedef __int8_t int8;
typedef __uint16_t uint16;
typedef int int32;
typedef unsigned int uint32;
typedef long int int64;
typedef unsigned long int uint64;

#include "Huffman.c"
#include "GlibIO.c"

#define USE_OMP
#define STRING_LEN (1024)

#define POW2_RANGE 30

int32 expbit, sigbit, expshift, SZnbit, SZnbyte;
uint32 expmask, sigmask;
int64 NtypeOfSample;
double _pow2_[POW2_RANGE*2+1];

#define MYEXP2(x) (_pow2_[((int)(x))+POW2_RANGE])

char inputfile[STRING_LEN], outputfile[STRING_LEN], logfile[STRING_LEN];

int8* Cbit;
uint8** C;
int32 MAX_COMPRESSION_ITR, * Nzero, * Cbyte;
uint32** sz_compressed_dat;
uint64** encoded_int_count;

float BREAK_RELATIVE_EPSILON, * trace_amp;
Trace* dat;

GlibHeader hdr;

FILE* flog;
char msgbuffer[STRING_LEN];
void printmsg(void){
    printf("%s", msgbuffer);
    fprintf(flog, "%s", msgbuffer);
    fflush(flog);
}

int64 index2d(int64 i, int64 j, int64 ni){ return ni*j+i; }

int64 lindex2d(int64* i, int64* j, int64 ni, int64 idx){
    *j = idx/ni;
    *i = idx-(*j)*ni;
    return idx;
}

void init_SZ_para(int32 _expbit, int32 _sigbit){
    int64 i, j;
    MYEXP2(0) = 1.0;
    for(i = 0; i<POW2_RANGE; i++){
        MYEXP2(i+1) = MYEXP2(i)*2.0;
        MYEXP2(-1-i) = MYEXP2(-i)*0.5;
    }
    expbit = _expbit;
    sigbit = _sigbit;
    expmask = (uint32)(MYEXP2(expbit)-1);
    sigmask = (uint32)(MYEXP2(sigbit)-1);
    SZnbit = expbit+sigbit+2;
    NtypeOfSample = (int64)MYEXP2(SZnbit);
    for(i = 0; i<POW2_RANGE; i++){
        j = (int64)MYEXP2(i);
        if(j>=SZnbit) break;
    }
    SZnbyte = (int32)(j/8);

    return;
}

float trace_maxabs(Trace t, int32 nt){
    int64 i, n;
    float maxv;
    maxv = 0.0;
    n = nt*18;
    for(i = 0; i<n; i++) if(maxv<fabsf(t[i])) maxv = fabsf(t[i]);
    return maxv;
}

int32 trace_nzero(Trace t, int32 nt){
    int32 nzero, it, j, flag;

    for(it = 0; it<nt; it++){
        flag = 0;
        for(j = 0; j<18; j++) if(t[index2d(it, j, nt)]!=0.0) flag |= 1;
        if(flag) break;
    }
    return it;
}

int32 calc_expshift(float maxv){
    double explimit;

    if(maxv==0.0)
        expshift = 0;
    else{
        explimit = floor(log2(fabs(maxv)))+1.0;
        expshift = (int32)(MYEXP2(expbit)-1.0-explimit);
    }
    return expshift;
}

float decode_float32(uint32 u){
    uint32 a, b, lim;
    int32 e;
    double s;

    lim = (uint32)MYEXP2(sigbit-1);
    a = u&sigmask;
    b = (u>>sigbit)&expmask;
    e = (int32)b-expshift;
    if(a<lim)
        s = ((double)a)*MYEXP2(e-sigbit+1);
    else
        s = ((double)a-MYEXP2(sigbit))*MYEXP2(e-sigbit+1);

    return (float)s;
}

uint32 encode_float32(float v, uint32 type){
    double _v, _fa;
    int32 e, sig;
    uint32 a, b, _u;

    _v = (double)v;
    if(fabs(_v)<decode_float32(1))
        e = -expshift;
    else
        e = floor(log2(fabs(_v)))+1.0;
    if(e<-expshift) e = -expshift;
    if((e+expshift+1)>=MYEXP2(expbit)) e = (MYEXP2(expbit)-1)-expshift;
    b = (uint32)(e+expshift);

    _fa = round(_v*MYEXP2(sigbit-e-1)); // v*2^(sigbit-e-1)
    if(_fa>=0.0){
        if(_fa>=MYEXP2(sigbit-1))
            _fa = MYEXP2(sigbit-1)-1.0;
    }
    else{
        _fa += MYEXP2(sigbit); // v * 2^(nsig-e-1) + 2^nsig
        if(_fa<MYEXP2(sigbit-1))
            _fa = MYEXP2(sigbit-1);
    }
    a = (uint32)_fa;
    _u = ((type&3)<<(expbit+sigbit))|((b&expmask)<<sigbit)|(a&sigmask);

    return _u;
}

float max_float(void){
    uint32 u, e, s;
    e = (uint32)(MYEXP2(expbit)-1.0);
    s = (uint32)(MYEXP2(sigbit-1)-1.0);
    u = ((e&expmask)<<sigbit)|(s&sigmask);
    return decode_float32(u);
}

double predict2p(double x1, double x2){ return 2*x2-x1; }

double predict3p(double x1, double x2, double x3){ return x1+3*(x3-x2); }

void trace_squeeze(uint32* encoded, float* decoded, float* max_abs_res,
    uint64* count, Trace t, float traceamp, int32 nzero, int32 nt){
    int64 i, j, idx, idx1, idx2, idx3;
    uint32 _type;
    float buf;
    double res, res1, res2, res3, res4, pre, pre1, pre2, pre3, pre4, _amp;

    _amp = (double)traceamp;
    for(i = 0; i<nt*18; i++){
        encoded[i] = 0;
        decoded[i] = 0.0;
    }
    *max_abs_res = 0.0;
    for(i = 0; i<NtypeOfSample; i++) count[i] = 0;
    for(j = 0; j<18; j++)
        for(i = nzero; i<nt; i++){
            idx = index2d(i, j, nt);
            idx1 = index2d(i-1, j, nt);
            idx2 = index2d(i-2, j, nt);
            idx3 = index2d(i-3, j, nt);
            switch(i){
                case 0:
                    pre1 = 0.0;
                    pre2 = 0.0;
                    pre3 = 0.0;
                    pre4 = 0.0;
                    break;
                case 1:
                    pre1 = decoded[idx1];
                    pre2 = predict2p(0.0, decoded[idx1]);
                    pre3 = predict3p(0.0, 0.0, decoded[idx1]);
                    pre4 = 0.0;
                    break;
                case 2:
                    pre1 = decoded[idx1];
                    pre2 = predict2p(decoded[idx2], decoded[idx1]);
                    pre3 = predict3p(0.0, decoded[idx2], decoded[idx1]);
                    pre4 = decoded[idx2];
                    break;
                default:
                    pre1 = decoded[idx1];
                    pre2 = predict2p(decoded[idx2], decoded[idx1]);
                    pre3 = predict3p(decoded[idx3], decoded[idx2], decoded[idx1]);
                    pre4 = decoded[idx2];
                    break;
            }
            res1 = (t[idx]-pre1)/_amp;
            res2 = (t[idx]-pre2)/_amp;
            res3 = (t[idx]-pre3)/_amp;
            res4 = (t[idx]-pre4)/_amp;
            if((fabs(res4)<fabs(res1))&&(fabs(res4)<fabs(res2))&&(fabs(res4)<fabs(res3))){
                res = res4;
                pre = pre4;
                _type = 3;
            }
            else if((fabs(res3)<fabs(res1))&&(fabs(res3)<fabs(res2))&&(fabs(res3)<fabs(res4))){
                res = res3;
                pre = pre3;
                _type = 2;
            }
            else if((fabs(res2)<fabs(res1))&&(fabs(res2)<fabs(res3))&&(fabs(res2)<fabs(res4))){
                res = res2;
                pre = pre2;
                _type = 1;
            }
            else{
                res = res1;
                pre = pre1;
                _type = 0;
            }
            if(fabs(res)>(*max_abs_res)) *max_abs_res = fabs(res);

            encoded[idx] = encode_float32(res, _type);
            count[encoded[idx]] += 1;
            decoded[idx] = pre+(double)decode_float32(encoded[idx])*_amp;
        }
    return;
}

void trace_huffman(uint8* compressed, HuffmanTree tree, uint32* hashtab,
    uint8** code, uint32* SZcode, int32* cbyte, int8* cbit, int32 nzero, int32 nt){
    uint32 ileaf;
    int32 it, icmp, icode;
    uint64 idx;
    *cbyte = -1;
    *cbit = 7;
    for(icmp = 0; icmp<18; icmp++)
        for(it = nzero; it<nt; it++){
            idx = index2d(it, icmp, nt);
            ileaf = hashtab[SZcode[idx]]-1;
            for(icode = 0; icode<tree.depth[ileaf]; icode++){
                *cbit += 1;
                if((*cbit)==8){
                    *cbyte += 1;
                    *cbit = 0;
                }
                if(code[ileaf][icode])
                    compressed[*cbyte] |= _BYTE_TRUES[*cbit];
                else
                    compressed[*cbyte] &= _BYTE_FALSES[*cbit];
            }
        }
    return;
}

void printhelp(void){
    printf("usage:\n");
    printf("\tcompress glibpath cglibpath [nexp] [nsig] [nthread] [maxitr] [breakeps]\n");
}

int main(int argc, char* argv[]){
    int8 ibuf8;
    uint8** codetable, uibuf8;
    uint16 uibuf16;
    int32 nt, nleaf, ibuf32, compress_itr, _nexp, _nsig, _nthread;
    uint32* hashtab, uibuf32;
    int64 i, j, k, idx, Ntrace;
    uint64* _weight_buffer, byte_written, uibuf64,
        part1pos, part2pos, part3pos, part4pos,
        hpos3, hpos2, * trace_start_pos;
    float maxval, lastmaxval, * maxval_array, ** decode_buffer, fbuf32;
    double fbuf64;
    HuffmanTree hftree;
    FILE* fout, * fp;

    if((argc>8)||(argc<3)){
        printhelp();
        exit(0);
    }
    if(argc<8)
        BREAK_RELATIVE_EPSILON = 1e-3;
    else
        BREAK_RELATIVE_EPSILON = atof(argv[7]);

    if(argc<7)
        MAX_COMPRESSION_ITR = 5;
    else
        MAX_COMPRESSION_ITR = atoi(argv[6]);

    if(argc<6)
        _nthread = 4;
    else
        _nthread = atoi(argv[5]);

    if(argc<5)
        _nsig = 10;
    else
        _nsig = atoi(argv[4]);

    if(argc<4)
        _nexp = 4;
    else
        _nexp = atoi(argv[3]);

    strcpy(inputfile, argv[1]);
    strcpy(outputfile, argv[2]);
    sprintf(logfile, "%s.log", inputfile);

    if(_nthread>1) omp_set_num_threads(_nthread);

    flog = fopen(logfile, "w");
    // init
    sprintf(msgbuffer, "init data\n"); printmsg();
    init_SZ_para(_nexp, _nsig);

    sprintf(msgbuffer, "  threads: %d, maxitr: %d, breakEps:%e\n",
        _nthread, MAX_COMPRESSION_ITR, BREAK_RELATIVE_EPSILON); printmsg();
    sprintf(msgbuffer, "  expbit: %d, expmask: 0x%04X\n", expbit, expmask); printmsg();
    sprintf(msgbuffer, "  sigbit: %d, sigmask: 0x%04X\n", sigbit, sigmask); printmsg();
    sprintf(msgbuffer, "  nbit: %d, nbyte: %d, ntype: %ld\n", SZnbit, SZnbyte, NtypeOfSample);
    printmsg();

    sprintf(msgbuffer, "read data\n"); printmsg();
    fp = fopen(inputfile, "r");
    glibio_read_head(&hdr, fp);
    glibio_read_data_all(&dat, hdr, fp);
    fclose(fp);
    Ntrace = ((int64)hdr.n[0])*((int64)hdr.n[1])*((int64)hdr.n[2]);
    nt = hdr.n[3];

    // allocate memory
    sprintf(msgbuffer, "allocate memory\n"); printmsg();
    Cbit = (int8*)calloc(Ntrace, sizeof(int8));
    C = (uint8**)calloc(Ntrace, sizeof(uint8*));
    for(i = 0; i<Ntrace; i++) C[i] = (uint8*)calloc(nt*18*sizeof(float), sizeof(uint8));
    Nzero = (int32*)calloc(Ntrace, sizeof(int32));
    Cbyte = (int32*)calloc(Ntrace, sizeof(int32));
    sz_compressed_dat = (uint32**)calloc(Ntrace, sizeof(uint32*));
    for(i = 0; i<Ntrace; i++) sz_compressed_dat[i] = (uint32*)calloc(nt*18, sizeof(uint32));
    encoded_int_count = (uint64**)calloc(Ntrace, sizeof(uint64*));
    for(i = 0; i<Ntrace; i++) encoded_int_count[i] = (uint64*)calloc(NtypeOfSample, sizeof(uint64));
    trace_amp = (float*)calloc(Ntrace, sizeof(float));

    _weight_buffer = (uint64*)calloc(NtypeOfSample, sizeof(uint64));
    maxval_array = (float*)calloc(Ntrace, sizeof(float));
    decode_buffer = (float**)calloc(Ntrace, sizeof(float*));
    for(i = 0; i<Ntrace; i++) decode_buffer[i] = (float*)calloc(nt*18, sizeof(float));


    // ! run length compress
#ifdef USE_OMP
#pragma omp parallel for private(i)
#endif
    for(i = 0; i<Ntrace; i++) Nzero[i] = trace_nzero(dat[i], nt);

#ifdef USE_OMP
#pragma omp parallel for private(i,j)
#endif
    for(i = 0; i<Ntrace; i++){
        trace_amp[i] = 0.0;
        for(j = 0; j<nt*18; j++) if(fabsf(dat[i][j])>trace_amp[i]) trace_amp[i] = fabsf(dat[i][j]);
    }
    maxval = 1.0;

    // ! SZ compress
    for(compress_itr = 0; compress_itr<MAX_COMPRESSION_ITR; compress_itr++){
        sprintf(msgbuffer, "%d step of SZ\n", compress_itr+1); printmsg();
        sprintf(msgbuffer, "  maxval: %f\n", maxval); printmsg();
        calc_expshift(maxval);
        sprintf(msgbuffer, "  expshift: %d\n", expshift); printmsg();

#ifdef USE_OMP
#pragma omp parallel for private(i)
#endif
        for(i = 0; i<Ntrace; i++)
            trace_squeeze(sz_compressed_dat[i], decode_buffer[i], &(maxval_array[i]),
            encoded_int_count[i], dat[i], trace_amp[i], Nzero[i], nt);
        lastmaxval = maxval;
        maxval = 0.0;
        for(i = 0; i<Ntrace; i++) if(maxval<maxval_array[i]) maxval = maxval_array[i];

        uibuf8 = 1;
        if(fabsf(maxval-lastmaxval)<(BREAK_RELATIVE_EPSILON*fabsf(lastmaxval)))
            uibuf8 &= 1;
        else
            uibuf8 &= 0;

        if(maxval>max_float())
            uibuf8 &= 0;
        else
            uibuf8 &= 1;

        if(uibuf8) break;
    }

    // ! Huffman Compress
    sprintf(msgbuffer, "prepare Huffman tree parameter\n"); printmsg();
    sprintf(msgbuffer, "  sum weight\n"); printmsg();
#ifdef USE_OMP
#pragma omp parallel for private(j)
#endif
    for(i = 0; i<NtypeOfSample; i++){
        _weight_buffer[i] = 0;
        for(j = 0; j<Ntrace; j++) _weight_buffer[i] += encoded_int_count[j][i];
    }
    sprintf(msgbuffer, "  count leaf\n"); printmsg();
    nleaf = 0;
    hashtab = (uint32*)calloc(NtypeOfSample, sizeof(uint32));
    for(i = 0; i<NtypeOfSample; i++)
        if(_weight_buffer[i]>0){
            hashtab[i] = nleaf+1;
            nleaf += 1;
        }
        else{
            hashtab[i] = 0;
        }
    sprintf(msgbuffer, "  nleaf: %d\n", nleaf); printmsg();
    sprintf(msgbuffer, "  allocate tree\n"); printmsg();
    huffman_alloc(&hftree, nleaf, sizeof(uint32));
    sprintf(msgbuffer, "  copy value\n"); printmsg();
    for(i = 0; i<NtypeOfSample; i++){
        j = hashtab[i];
        if(j>0){
            ((uint32*)hftree.value)[j-1] = i;
            hftree.weight[j-1] = (int64)_weight_buffer[i];
        }
    }
    sprintf(msgbuffer, "  hftree.Nleaf: %d, hftree.Nnode: %d\n", hftree.Nleaf, hftree.Nnode);
    printmsg();
    sprintf(msgbuffer, "construct huffman tree\n"); printmsg();
    sprintf(msgbuffer, "  before:\n"); printmsg();
    huffman_fprintf(flog, hftree);

    sprintf(msgbuffer, "  construct:\n"); printmsg();
    huffman_construct(&hftree);

    sprintf(msgbuffer, "  after:\n"); printmsg();
    huffman_fprintf(flog, hftree);

    sprintf(msgbuffer, "make encode table\n"); printmsg();
    codetable = (uint8**)calloc(hftree.Nleaf, sizeof(uint8*));
    huffman_encode(hftree, codetable);
    sprintf(msgbuffer, "  hftree.Nleaf: %d, hftree.Nnode: %d\n", hftree.Nleaf, hftree.Nnode);
    printmsg();

    sprintf(msgbuffer, "huffman encode\n"); printmsg();

#ifdef USE_OMP
#pragma omp parallel for private(i)
#endif
    for(i = 0; i<Ntrace; i++)
        trace_huffman(C[i], hftree, hashtab, codetable, sz_compressed_dat[i],
        &(Cbyte[i]), &(Cbit[i]), Nzero[i], nt);

    // ! write
    sprintf(msgbuffer, "write to file\n"); printmsg();
    fout = fopen(outputfile, "w"); fseek(fout, 0, SEEK_SET);
    byte_written = 0;
    uibuf8 = 0x71; byte_written += fwrite(&uibuf8, sizeof(uint8), 1, fout)*sizeof(uint8);
    uibuf64 = 0;
    for(i = 0; i<4; i++)
        byte_written += fwrite(&uibuf64, sizeof(uint64), 1, fout)*
        sizeof(uint64);
    // srclat srclon, srcdep
    fbuf32 = 0; for(i = 0; i<3; i++)
        byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float);
    // nt
    byte_written += fwrite(&nt, sizeof(int32), 1, fout)*sizeof(int32);
    // dt
    fbuf32 = hdr.t[1]-hdr.t[0];
    byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float);
    // stf
    fbuf32 = 0.0;
    for(i = 0; i<nt; i++)
        byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float);
    fflush(fout);
    part1pos = byte_written;
    // nx, ny, nz
    byte_written += fwrite(hdr.n, sizeof(int32), 3, fout)*sizeof(int32);
    byte_written += fwrite(&(hdr.x[0]), sizeof(float), 1, fout)*sizeof(float); // x0
    byte_written += fwrite(&(hdr.y[0]), sizeof(float), 1, fout)*sizeof(float); // y0
    byte_written += fwrite(&(hdr.z[0]), sizeof(float), 1, fout)*sizeof(float); // z0
    fbuf32 = hdr.x[1]-hdr.x[0]; byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float); // dx
    fbuf32 = hdr.y[1]-hdr.y[0]; byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float); // dy
    fbuf32 = hdr.z[1]-hdr.z[0]; byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float); // dz
    hpos3 = byte_written;
    uibuf64 = 0;
    byte_written += fwrite(&uibuf64, sizeof(uint64), 1, fout)*sizeof(uint64);
    // air flag
    uibuf8 = 0xff;
    for(i = 0; i<floor((Ntrace-1)/8.0)+1; i++)
        byte_written += fwrite(&uibuf8, sizeof(uint8), 1, fout)*sizeof(uint8);
    fflush(fout);
    hpos2 = byte_written;
    // trace start position
    trace_start_pos = (uint64*)calloc(Ntrace, sizeof(uint64));
    byte_written += fwrite(trace_start_pos, sizeof(uint64), Ntrace, fout)*sizeof(uint64);
    part2pos = byte_written;
    ibuf8 = expbit; byte_written += fwrite(&ibuf8, sizeof(int8), 1, fout)*sizeof(int8);
    ibuf8 = sigbit; byte_written += fwrite(&ibuf8, sizeof(int8), 1, fout)*sizeof(int8);
    ibuf32 = expshift; byte_written += fwrite(&ibuf32, sizeof(int32), 1, fout)*sizeof(int32);
    fbuf64 = maxval; byte_written += fwrite(&fbuf64, sizeof(double), 1, fout)*sizeof(double);
    byte_written += fwrite(&(hftree.Nleaf), sizeof(int32), 1, fout)*sizeof(int32);
    for(i = 0; i<hftree.Nleaf; i++){
        switch(SZnbyte){
            case 1:
                uibuf8 = ((uint32*)hftree.value)[i];
                byte_written += fwrite(&uibuf8, sizeof(uint8), 1, fout)*sizeof(uint8);
                break;
            case 2:
                uibuf16 = ((uint32*)hftree.value)[i];
                byte_written += fwrite(&uibuf16, sizeof(uint16), 1, fout)*sizeof(uint16);
                break;
            case 4:
                uibuf32 = ((uint32*)hftree.value)[i];
                byte_written += fwrite(&uibuf32, sizeof(uint32), 1, fout)*sizeof(uint32);
                break;
            case 8:
                uibuf64 = ((uint32*)hftree.value)[i];
                byte_written += fwrite(&uibuf64, sizeof(uint64), 1, fout)*sizeof(uint64);
                break;
        }
    }
    byte_written += fwrite(&(hftree.Nnode), sizeof(int32), 1, fout)*sizeof(int32);
    for(i = 0; i<hftree.Nnode; i++){
        ibuf32 = hftree.left[i];
        byte_written += fwrite(&ibuf32, sizeof(int32), 1, fout)*sizeof(int32);
    }
    for(i = 0; i<hftree.Nnode; i++){
        ibuf32 = hftree.right[i];
        byte_written += fwrite(&ibuf32, sizeof(int32), 1, fout)*sizeof(int32);
    }
    fflush(fout);
    part3pos = byte_written;
    uibuf8 = 0; for(i = 0; i<1024; i++) byte_written += fwrite(&uibuf8, sizeof(uint8), 1, fout)*sizeof(uint8);
    fflush(fout);
    part4pos = byte_written;
    for(idx = 0; idx<Ntrace; idx++){
        trace_start_pos[idx] = byte_written;
        fbuf32 = 0.0;
        // tp
        byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float);
        // ts
        byte_written += fwrite(&fbuf32, sizeof(float), 1, fout)*sizeof(float);
        // nzero
        byte_written += fwrite(&(Nzero[idx]), sizeof(int32), 1, fout)*sizeof(int32);
        // amp
        byte_written += fwrite(&(trace_amp[idx]), sizeof(float), 1, fout)*sizeof(float);
        // nbyte
        ibuf32 = Cbyte[idx]+1; byte_written += fwrite(&ibuf32, sizeof(int32), 1, fout)*sizeof(int32);
        // nbit
        ibuf8 = Cbit[idx]+1; byte_written += fwrite(&ibuf8, sizeof(int8), 1, fout)*sizeof(int8);
        // data
        byte_written += fwrite(C[idx], sizeof(uint8), Cbyte[idx]+1, fout)*sizeof(uint8);
        fflush(fout);
    }
    fseek(fout, 1, SEEK_SET);
    fwrite(&part1pos, sizeof(uint64), 1, fout);
    fwrite(&part2pos, sizeof(uint64), 1, fout);
    fwrite(&part3pos, sizeof(uint64), 1, fout);
    fwrite(&part4pos, sizeof(uint64), 1, fout);

    fflush(fout);
    fseek(fout, hpos3, SEEK_SET);
    fwrite(&hpos2, sizeof(uint64), 1, fout);

    fflush(fout);
    fseek(fout, hpos2, SEEK_SET);
    fwrite(trace_start_pos, sizeof(uint64), Ntrace, fout);
    fflush(fout);

    // free memory
    sprintf(msgbuffer, "free memory\n"); printmsg();
    free(Cbit);
    for(i = 0; i<Ntrace; i++) free(C[i]); free(C);
    free(Nzero);
    free(Cbyte);
    for(i = 0; i<Ntrace; i++) free(sz_compressed_dat[i]); free(sz_compressed_dat);
    for(i = 0; i<Ntrace; i++) free(encoded_int_count[i]); free(encoded_int_count);
    free(trace_amp);
    glibio_free_data(&dat, hdr);
    glibio_free_head(&hdr);
    for(i = 0; i<hftree.Nleaf; i++) free(codetable[i]);
    free(codetable);
    free(hashtab);
    free(_weight_buffer);
    free(trace_start_pos);
    free(maxval_array);
    for(i = 0; i<Ntrace; i++) free(decode_buffer[i]); free(decode_buffer);
    huffman_free(&hftree);
    fclose(flog);
    fclose(fout);

    return 0;
}
