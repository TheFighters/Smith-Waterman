/***********************************************************************
 * Smith–Waterman algorithm
 * Purpose:     Local alignment of nucleotide or protein sequences
 * Authors:     Daniel Holanda, Hanoch Griner, Taynara Pinheiro
 ***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

/*--------------------------------------------------------------------
 * Text Tweaks
 */
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */
/* End of text tweaks */

/*--------------------------------------------------------------------
 * Constants
 */
#define UP 1
#define LEFT 2
#define DIAGONAL 3
#define PATH -1

/* End of constants */


/*--------------------------------------------------------------------
 * Functions Prototypes
 */
void similarityScore(int i, int j, int* H, int* P, int* B, int* maxPos);
int matchMissmatchScore(int i, int j);
void backtrack(int* B, int maxPos);
void printMatrix(int* matrix);
void printPredecessorMatrix(int* matrix, int* B);
void printPositionMatrix(void);
/* End of prototypes */


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

    //Allocates predecessor matrix B
    int *B;
    B = malloc(m * n * sizeof(int));

    //Start position for backtrack
    int maxPos = 0;

    //Calculates the similarity matrix
    int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) { //Columns
            similarityScore(i, j, H, P, B, &maxPos);
        }
    }

    backtrack(B, maxPos);

    printf("\nSimilarity Matrix:\n");
    printMatrix(H);

    printf("\nPredecessor Matrix:\n");
    printPredecessorMatrix(P, B);
    
    printf("\nMaxPos: %d\n", maxPos);

    printf("\nPosition Matrix:\n");
    printPositionMatrix();

    printf("\nPredecessor Backtrack Matrix:\n");
    printMatrix(B);
    //Frees similarity matrixs
    free(H);
    free(P);
    free(B);

    return 0;
}  /* End of main */


/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate  the maximum Similarity-Score H(i,j)
 */
void similarityScore(int i, int j, int* H, int* P, int* B, int* maxPos) {

    int up, left, diag;

    //Stores index of element
    int index = m * i + j;

    //Get element above
    if (i == 0) {
        up = gapScore;
    }
    else {
        up = H[index - m] + gapScore;
    }

    //Get element on the left
    if (j == 0) {
        left = gapScore;
    }
    else {
        left = H[index - 1] + gapScore;
    }

    //Get element on the diagonal
    if (j == 0 || i == 0) {
        diag = matchMissmatchScore(i, j);
    }
    else {
        diag = H[index - m - 1] + matchMissmatchScore(i, j);
    }

    //Calculates the maximum
    int max = 0;
    int pred, predPos;
    if (up > max) {
        max = up;
        pred = UP;
        predPos = index - m;
    }
    if (left > max) {
        max = left;
        pred = LEFT;
        predPos = index - 1;
    }
    if (diag > max) {
        max = diag;
        pred = DIAGONAL;
        predPos = index - m - 1;
    }

    //Inserts the value in the similarity and predecessor matrixes
    H[index] = max;
    P[index] = pred;

    //Check if predecessor position is valid
    //first line has to go left and first col has to go up 
    if ((i != 0 && j != 0) || i == 0 && pred == LEFT || j == 0 && pred == UP)
        B[index] = predPos;
    else
        B[index] = 0;

    //Updates maximum score to be used as seed on backtrack 
    if (max > H[*maxPos]) {
        *maxPos = index;
    }

}  /* End of similarityScore */


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
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(int* B, int maxPos) {
    //hold maxPos value
    int predPos;

    //backtrack from maxPos to startPos = 0 
    do {
        predPos = B[maxPos];
        B[maxPos] = PATH;
        maxPos = predPos;
    } while(B[maxPos] != 0);
     
    //sets start '0' to '-1'
    B[maxPos] = PATH;

}

/*--------------------------------------------------------------------
 * Function:    printMatrix
 * Purpose:     Print Matrix
 */
void printMatrix(int* matrix) {
    int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            printf("%d\t", matrix[m * i + j]);
        }
        printf("\n");
    }

}  /* End of printMatrix */

/*--------------------------------------------------------------------
 * Function:    printPredecessorMatrix
 * Purpose:     Print predecessor matrix
 */
void printPredecessorMatrix(int* matrix, int* B) {
    int i, j, index, path;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            index = m * i + j;
            path = B[index] + 1;
            if (matrix[index] == UP) {
                if (path) printf("↑ ");
                else printf(BOLDGREEN "↑ " RESET);
            }
            else if (matrix[index] == LEFT) {
                if (path) printf("← ");
                else printf(BOLDGREEN "← " RESET);
            }
            else if (matrix[index] == DIAGONAL) {
                if (path) printf("↖ ");
                else printf(BOLDGREEN "↖ " RESET);
            }
            else {
                printf("- ");
            }
        }
        printf("\n");
    }

}  /* End of printPredecessorMatrix */

/*--------------------------------------------------------------------
 * Function:    printPositionMatrix
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
void printPositionMatrix(void) {
    int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            printf("%d\t", m*i+j);
        }
        printf("\n");
    }
}  /* End of printPositionMatrix */



/*--------------------------------------------------------------------
 * External References:
 * http://vlab.amrita.edu/?sub=3&brch=274&sim=1433&cnt=1
 */
