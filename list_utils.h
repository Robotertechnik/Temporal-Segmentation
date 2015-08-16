#ifndef LIST_UTILS
#define LIST_UTILS


#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>

typedef struct {
  int Max;          /* taille max de la Lifo */
  int Sp;           /* index de pile (pointe la 1ere case libre) */
  int Pts[1];
} Lifo;


Lifo * CreeLifoVide(int taillemax) {
  Lifo * L = (Lifo *)calloc(1,sizeof(Lifo) + sizeof(int32_t) * (taillemax-1));
 if (L == NULL) {   
    fprintf(stderr, "CreeLifoVide() : malloc failed : %ld bytes\n", 
            sizeof(Lifo) + sizeof(int32_t) * (taillemax-1));
    return NULL;
  }
  L->Max = taillemax;
  L->Sp = 0;
  return L;
}

int LifoVide(Lifo * L) {
  return (L->Sp == 0);
}

int LifoPop(Lifo * L) {
  if (L->Sp == 0) {
    fprintf(stderr, "erreur Lifo vide\n");
    exit(1);
  }
  L->Sp -= 1;
  return L->Pts[L->Sp];
}

void LifoPush(Lifo * L, int V) {
  if (L->Sp > L->Max - 1) {
    fprintf(stderr, "erreur Lifo pleine\n");
    exit(1);
  }
  L->Pts[L->Sp] = V;
  L->Sp += 1;
}

void LifoPrint(Lifo * L) {
  int i;
  if (LifoVide(L)) {printf("[]"); return;}
  printf("[ ");
  for (i = 0; i < L->Sp; i++)
    printf("%ld ", L->Pts[i]);
  printf("]");
}

void LifoTermine(Lifo * L) {
  free(L);
}


#endif
