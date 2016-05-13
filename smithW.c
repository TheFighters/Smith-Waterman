/***********************************************************************
 * Smith–Waterman algorithm
 * Purpose:     Local alignment of nucleotide or protein sequences
 * Authors:     Daniel Holanda, Hanoch Griner, Taynara Pinheiro
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void similarityScore(int i, int j, int* H, int* P);
void printSimilarityMatrix(int* matrix);
void printPredecessorMatrix(int* matrix);
int matchMissmatchScore(int i, int j);


/*--------------------------------------------------------------------
 * Global Variables
 */
//Defines size of strings to be compared
int m = 11; //Columns - Size of string a
int n = 7;  //Lines - Size of string b

//Defines scores
int matchScore = 5;
int missmatchScore = -3;
int gapScore = -4;

//Strings over the Alphabet Sigma
char a[] = {'C', 'G', 'T', 'G', 'A', 'A', 'T', 'T', 'C', 'A', 'T'};
char b[] = {'G', 'A', 'C', 'T', 'T', 'A', 'C'};

/* End of global variables */

/*--------------------------------------------------------------------
 * Function:    main
 */
int main(int argc, char* argv[]) {

    //Allocates similarity matrix H
    int *H;
    H = malloc(m * n * sizeof(int));

    //Allocates predecessor matrix P
    int *P;
    P = malloc(m * n * sizeof(int));

    //Calculates the similarity matrix
    int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) { //Columns
            similarityScore(i, j, H, P);
        }
    }

    printf("\nSimilarity Matrix:\n");
    printSimilarityMatrix(H);
    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P);

    //Frees similarity matrix H
    free(H);

    return 0;
}  /* End of main */


/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate  the maximum Similarity-Score H(i,j)
 */
void similarityScore(int i, int j, int* H, int* P) {

    int up, left, diag;

    //Get element above
    if (i == 0) {
        up = gapScore;
    }
    else {
        up = H[m * (i - 1) + j] + gapScore;
    }

    //Get element on the left
    if (j == 0) {
        left = gapScore;
    }
    else {
        left = H[m * i + j - 1] + gapScore;
    }

    //Get element on the diagonal
    if (j == 0 || i == 0) {
        diag = matchMissmatchScore(i, j);
    }
    else {
        diag = H[m * (i - 1) + j - 1] + matchMissmatchScore(i, j);
    }

    //Calculates the maximum
    int max = 0;
    int pred = 0;
    if (up > max) {
        max = up;
        pred = 1;
    }
    if (left > max) {
        max = left;
        pred = 2;
    }
    if (diag > max) {
        max = diag;
        pred = 3;
    }

    //Inserts the value in the similarity adn predecessor matrixes
    H[m * i + j] = max;
    P[m * i + j] = pred;

}  /* End of printSimilarityMatrix */

/*--------------------------------------------------------------------
 * Function:    printSimilarityMatrix
 * Purpose:     Print similarity matrix
 */
void printSimilarityMatrix(int* matrix) {
    int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            printf("%d\t", matrix[m * i + j]);
        }
        printf("\n");
    }

}  /* End of similarityScore */

/*--------------------------------------------------------------------
 * Function:    printPredecessorMatrix
 * Purpose:     Print predecessor matrix
 */
void printPredecessorMatrix(int* matrix) {
    int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            if (matrix[m * i + j] == 1) {
                printf("↑ ");
            }
            else if (matrix[m * i + j] == 2) {
                printf("← ");
            }
            else if (matrix[m * i + j] == 3) {
                printf("↖ ");
            }
            else {
                printf("- ");
            }
        }
        printf("\n");
    }

}  /* End of printPredecessorMatrix */

/*--------------------------------------------------------------------
 * Function:    matchMissmatchScore
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
int matchMissmatchScore(int i, int j) {
    if (a[j] == b[i])
        return matchScore;
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */

/*--------------------------------------------------------------------
 * External References:
 * http://vlab.amrita.edu/?sub=3&brch=274&sim=1433&cnt=1
 */
