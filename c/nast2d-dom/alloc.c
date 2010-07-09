#include <stdlib.h>

/* Allocate memory for a rows*cols array of doubles.
 * The elements within a column are contiguous in memory, and columns
 * themselves are also contiguous in memory.
 */
double **alloc_doublematrix(int cols, int rows)
{
    int i;
    double **m;
    if ((m = (double**) malloc(cols*sizeof(double*))) == NULL) {
        return NULL;
    }
    double *els = (double *) calloc(rows*cols, sizeof(double));
    if (els == NULL) {
        return NULL;
    } 
    for (i = 0; i < cols; i++) {
        m[i] = &els[rows * i];
    }
    return m;
} 

/* Allocate memory for a rows*cols array of chars. */
char **alloc_charmatrix(int cols, int rows)
{
    int i;
    char **m;
    if ((m = (char**) malloc(cols*sizeof(char*))) == NULL) {
        return NULL;
    }
    char *els = (char *) malloc(rows*cols*sizeof(char));
    if (els == NULL) {
        return NULL;
    } 
    for (i = 0; i < cols; i++) {
        m[i] = &els[rows * i];
    }
    return m;
} 

/* Free the memory of a matrix allocated with alloc_{double|char}matrix*/
void free_matrix(void *m)
{
    void **els = (void **) m;
    free(els[0]);   /* Deallocate the block of array elements */
    free(m);        /* Deallocate the block of column pointers */
}
