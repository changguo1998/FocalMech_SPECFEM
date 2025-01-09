#pragma GCC diagnostic ignored "-Wunused-result"

#ifndef __GLIBIO_C__
#define __GLIBIO_C__

#include <stdio.h>
#include <stdlib.h>

typedef int int32;
typedef unsigned int uint32;
typedef long int int64;
typedef unsigned long int uint64;

typedef struct _GlibHeader{
    float rt;
    int32 n[4];
    float* x, * y, * z, * t;
} GlibHeader;

typedef float* Trace;

int32 glibio_read_head(GlibHeader* hdr, FILE* fp){
    fseek(fp, 0, SEEK_SET);
    fread(&(hdr->rt), sizeof(float), 1, fp);
    fread(hdr->n, sizeof(int32), 4, fp);
    hdr->x = (float*)malloc(sizeof(float)*hdr->n[0]);
    hdr->y = (float*)malloc(sizeof(float)*hdr->n[1]);
    hdr->z = (float*)malloc(sizeof(float)*hdr->n[2]);
    hdr->t = (float*)malloc(sizeof(float)*hdr->n[3]);
    fread(hdr->x, sizeof(float), hdr->n[0], fp);
    fread(hdr->y, sizeof(float), hdr->n[1], fp);
    fread(hdr->z, sizeof(float), hdr->n[2], fp);
    fread(hdr->t, sizeof(float), hdr->n[3], fp);
    return 0;
}

void glibio_write_head(GlibHeader hdr, FILE* _fp){
    fwrite(&(hdr.rt), sizeof(float), 1, _fp);
    fwrite(hdr.n, sizeof(int32), 4, _fp);
    fwrite(hdr.x, sizeof(float), hdr.n[0], _fp);
    fwrite(hdr.y, sizeof(float), hdr.n[1], _fp);
    fwrite(hdr.z, sizeof(float), hdr.n[2], _fp);
    fwrite(hdr.t, sizeof(float), hdr.n[3], _fp);
    return;
}

void glibio_free_head(GlibHeader* hdr){
    free(hdr->x); hdr->x = NULL;
    free(hdr->y); hdr->y = NULL;
    free(hdr->z); hdr->z = NULL;
    free(hdr->t); hdr->t = NULL;
    return;
}

void glibio_read_trace(Trace* tr, int32 nt, FILE* fp){
    *tr = (Trace)malloc(sizeof(float)*nt*18);
    fread(*tr, sizeof(float), nt*18, fp);
    return;
}

void glibio_write_trace(Trace tr, int32 nt, FILE* fp){
    fwrite(tr, sizeof(float), nt*18, fp);
    return;
}

void glibio_free_trace(Trace* tr){
    free(*tr); *tr = NULL;
    return;
}

void glibio_read_data_all(Trace** buf, GlibHeader hdr, FILE* fp){
    int64 i, Ntrace;
    Ntrace = hdr.n[0]*hdr.n[1]*hdr.n[2];
    *buf = (Trace*)malloc(sizeof(Trace*)*Ntrace);
    for(i = 0; i<Ntrace; i++) glibio_read_trace(&((*buf)[i]), hdr.n[3], fp);
}

void glibio_write_data_all(Trace* dat, GlibHeader hdr, FILE* fp){
    int64 i, Ntrace;
    Ntrace = hdr.n[0]*hdr.n[1]*hdr.n[2];
    for(i = 0; i<Ntrace; i++) glibio_write_trace(dat[i], hdr.n[3], fp);
    return;
}

void glibio_free_data(Trace** buf, GlibHeader hdr){
    int64 i, Ntrace;
    Ntrace = hdr.n[0]*hdr.n[1]*hdr.n[2];
    for(i = 0; i<Ntrace; i++) glibio_free_trace(&((*buf)[i]));
    free(*buf);
    *buf = NULL;
    return;
}

#endif