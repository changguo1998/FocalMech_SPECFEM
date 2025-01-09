#ifndef __HUFFMAN__
#define __HUFFMAN__

// #define __HUFFMAN_DEBUG__

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>

typedef int int32;
typedef unsigned int uint32;
typedef long int int64;
typedef unsigned long int uint64;

__uint8_t _BYTE_TRUES[8] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 },
_BYTE_FALSES[8] = { 0x7F, 0xBF, 0xDF, 0xEF, 0xF7, 0xFB, 0xFD, 0xFE };

/*
 * QueueNode
 *     - struct QueueNode* next;
 *     - int64 v;
*/
typedef struct _queue_node{
    struct _queue_node* next;
    int64 v;
} QueueNode;

typedef struct _queue{
    QueueNode* head;
    QueueNode* tail;
    QueueNode* list;
} Queue;

void queue_init(Queue* q){
    q->head = NULL;
    q->tail = NULL;
    q->list = NULL;
    return;
}

void queue_push(Queue* q, int64 v){
    QueueNode* new_node;
    new_node = (QueueNode*)malloc(sizeof(QueueNode));
    new_node->next = NULL;
    new_node->v = v;
    if(q->list==NULL){
        q->list = new_node;
        q->tail = new_node;
        q->head = new_node;
    }
    else{
        (q->tail)->next = new_node;
        q->tail = new_node;
    }
    return;
}

int64 queue_pop(Queue* q, int64* _v){
    QueueNode* _h;
    if(q->list==NULL) return 1;

    _h = q->head;
    *_v = _h->v;
    q->head = _h->next;
    q->list = _h->next;
    free(_h);
    return 0;
}

int64 queue_popindex(Queue* q, int64* _v, int64 _idx){
    int64 i;
    QueueNode* _p;
    QueueNode* _q;
    if(_idx==0) return queue_pop(q, _v);

    _q = q->head;
    for(i = 0; (i<_idx)&&(_q!=NULL); i++){
        _p = _q;
        _q = _q->next;
    }

    if(i!=_idx) return 1;

    *_v = _q->v;
    _p->next = _q->next;
    free(_q);

    return 0;
}

int64 queue_popvalue(Queue* q, int64 _v){
    int64 i;
    QueueNode* _p;
    QueueNode* _q;
    if((q->head)->v==_v) return queue_pop(q, &i);
    _q = q->head;
    while((_q->v!=_v)&&(_q!=NULL)){
        _p = _q;
        _q = _q->next;
    }
    if(_q==NULL) return 1;
    if(_q->v!=_v) return 1;
    if(_q->next==NULL) q->tail = _p;
    _p->next = _q->next;
    free(_q);

    return 0;
}

typedef struct _HuffmanTree{
    int32 Nleaf;
    int32 Nnode;
    void* value;
    uint64* weight;
    int32* left;
    int32* right;
    int32* parent;
    int32* depth;
    int32* isright;
} HuffmanTree;

void huffman_init(HuffmanTree* t){
    t->Nleaf = 0;
    t->Nnode = 0;
    t->value = NULL;
    t->weight = NULL;
    t->left = NULL;
    t->right = NULL;
    t->parent = NULL;
    t->depth = NULL;
    t->isright = NULL;
    return;
}

void huffman_alloc(HuffmanTree* t, int32 Nleaf, size_t value_byte){
    int32 Nnode;
    Nnode = 2*Nleaf-1;
    t->Nleaf = Nleaf;
    t->Nnode = Nnode;
    t->value = calloc(Nnode, value_byte);
    t->weight = (uint64*)calloc(Nnode, sizeof(uint64));
    t->left = (int32*)calloc(Nnode, sizeof(int32));
    t->right = (int32*)calloc(Nnode, sizeof(int32));
    t->parent = (int32*)calloc(Nnode, sizeof(int32));
    t->depth = (int32*)calloc(Nnode, sizeof(int32));
    t->isright = (int32*)calloc(Nnode, sizeof(int32));
    return;
}

void huffman_construct(HuffmanTree* t){
    Queue nodelist;
    QueueNode* _qn;
    int64 i, inew, i1, i2, q1, q2;
    uint64  w1, w2, w0;

#ifdef __HUFFMAN_DEBUG__
    printf("init nodelist\n");
#endif
    queue_init(&nodelist);
    // construct init list
#ifdef __HUFFMAN_DEBUG__
    printf("push nodelist\n");
#endif
    for(i = 0; i<t->Nleaf; i++) queue_push(&nodelist, i+1);
#ifdef __HUFFMAN_DEBUG__
    _qn = nodelist.head;
    while(_qn!=NULL){
        printf("%d ", _qn->v);
        _qn = _qn->next;
    }
    printf("\n");
#endif
    inew = t->Nleaf+1;
    while(nodelist.head->next!=NULL){
        // find max
        _qn = nodelist.head;
        w0 = 0;
        while(_qn!=NULL){
            if(w0<t->weight[_qn->v-1]){
                w0 = t->weight[_qn->v-1];
            }
            _qn = _qn->next;
        }
        // find min 1
        _qn = nodelist.head;
        w1 = w0+1;
        i = 0;
        i1 = 0;
        q1 = 0;
        do{
            if(w1>t->weight[_qn->v-1]){
                q1 = i;
                i1 = _qn->v;
                w1 = t->weight[_qn->v-1];
            }
            i++;
            _qn = _qn->next;
        } while(_qn!=NULL);
        // find min 2
        _qn = nodelist.head;
        w2 = w0+1;
        i = 0;
        i2 = 0;
        q2 = 0;
        do{
            if((w2>t->weight[_qn->v-1])&&(i!=q1)){
                q2 = i;
                i2 = _qn->v;
                w2 = t->weight[_qn->v-1];
            }
            i++;
            _qn = _qn->next;
        } while(_qn!=NULL);
#ifdef __HUFFMAN_DEBUG__
        printf("inew: %d, i1: %d, i2: %d\n", inew, i1, i2);
#endif
#ifdef __HUFFMAN_DEBUG__
        _qn = nodelist.head;
        while(_qn!=NULL){
            printf("%d ", _qn->v);
            _qn = _qn->next;
        }
        printf("\n");
#endif
        if(queue_popvalue(&nodelist, i1)>0) exit(1);
#ifdef __HUFFMAN_DEBUG__
        _qn = nodelist.head;
        while(_qn!=NULL){
            printf("%d ", _qn->v);
            _qn = _qn->next;
        }
        printf("\n");
#endif
        if(queue_popvalue(&nodelist, i2)>0) exit(1);
#ifdef __HUFFMAN_DEBUG__
        _qn = nodelist.head;
        while(_qn!=NULL){
            printf("%d ", _qn->v);
            _qn = _qn->next;
        }
        printf("\n");
#endif
        t->parent[i1-1] = inew;
        t->isright[i1-1] = 0;
        t->parent[i2-1] = inew;
        t->isright[i2-1] = 1;
        t->left[inew-1] = i1;
        t->right[inew-1] = i2;
        t->weight[inew-1] = t->weight[i1-1]+t->weight[i2-1];
        queue_push(&nodelist, inew);
#ifdef __HUFFMAN_DEBUG__
        _qn = nodelist.head;
        while(_qn!=NULL){
            printf("%d ", _qn->v);
            _qn = _qn->next;
        }
        printf("\n");
#endif
        inew++;
    }

    // clean nodelist
    while(nodelist.head!=NULL) queue_pop(&nodelist, &i);
    // calculate depth
    queue_push(&nodelist, t->Nnode);
    t->depth[t->Nnode-1] = 0;
    while(nodelist.list!=NULL){
        queue_pop(&nodelist, &inew);
        w0 = t->depth[inew-1];
        i1 = t->left[inew-1];
        if(i1>0){
            queue_push(&nodelist, i1);
            t->depth[i1-1] = w0+1;
        }
        i2 = t->right[inew-1];
        if(i2>0){
            queue_push(&nodelist, i2);
            t->depth[i2-1] = w0+1;
        }
    }

    return;
}

void huffman_encode(HuffmanTree t, __uint8_t** c){
    int64 maxdep, i, j, inode;
    maxdep = 0;
    for(i = 0; i<t.Nleaf; i++) if(maxdep<t.depth[i]) maxdep = t.depth[i];
#ifdef __HUFFMAN_DEBUG__
    printf("maxdep: %d\n", maxdep);
#endif
    // *c = (__uint8_t**)calloc(maxdep*t.Nleaf, sizeof(__uint8_t));
    for(i = 0; i<t.Nleaf; i++)
        c[i] = (__uint8_t*)calloc(maxdep, sizeof(__uint8_t));
    // (*c)[0][0] = 0x00;
#ifdef __HUFFMAN_DEBUG__
    printf("allocated\n");
#endif
    for(i = 0; i<t.Nleaf; i++){
#ifdef __HUFFMAN_DEBUG__
        printf("=i= %d\n", i);
#endif
        inode = i+1;
        j = t.depth[i]-1;
        while(inode!=t.Nnode){
#ifdef __HUFFMAN_DEBUG__
            printf("%d ", inode);
#endif
            c[i][j] = (__uint8_t)(t.isright[inode-1]);
            inode = t.parent[inode-1];
            j -= 1;
        }
#ifdef __HUFFMAN_DEBUG__
        printf("\n");
#endif
    }
    return;
}

void huffman_printf(HuffmanTree t){
    int64 i;
    for(i = 0; i<t.Nnode; i++)
        printf("%d %ld %d %d %d %d %d\n", ((uint32*)t.value)[i], t.weight[i], \
        t.parent[i], t.left[i], t.right[i], t.isright[i], t.depth[i]);
    return;
}

void huffman_fprintf(FILE* fp, HuffmanTree t){
    int64 i;
    for(i = 0; i<t.Nnode; i++)
        fprintf(fp, "%d %ld %d %d %d %d %d\n", ((uint32*)t.value)[i], t.weight[i], \
        t.parent[i], t.left[i], t.right[i], t.isright[i], t.depth[i]);
    return;
}

void huffman_free(HuffmanTree* t){
    free(t->value); t->value = NULL;
    free(t->weight); t->weight = NULL;
    free(t->left); t->left = NULL;
    free(t->right); t->right = NULL;
    free(t->parent); t->parent = NULL;
    free(t->depth); t->depth = NULL;
    free(t->isright); t->isright = NULL;
    return;
}

#endif
