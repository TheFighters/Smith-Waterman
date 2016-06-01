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
#define NONE 0
#define UP 1
#define LEFT 2
#define DIAGONAL 3
#define PATH -1
#define asd 12


//Defines size of strings to be compared
#define m 10000        //Columns - Size of string a
#define n 10000           //Lines - Size of string b

//Defines scores
#define matchScore 5;
#define missmatchScore -3;
#define gapScore -4;

/*--------------------------------------------------------------------
* Helpers
*/
#define min(x, y) (((x) < (y)) ? (x) : (y))
#define max(a,b) ((a) > (b) ? a : b)

/* End of Helpers */

/*--------------------------------------------------------------------
 * Functions Prototypes
 */
void similarityScore(long int i, long int j, int* H, int* P, int* B, long int* maxPos);
long int matchMissmatchScore(long int i, long int j);
void backtrack(int* B, long int maxPos);
void printMatrix(int* matrix);
void printPredecessorMatrix(int* matrix, int* B);
void printPositionMatrix(void);
long int nElement(long int i);
void calcIndex(long int i, long int j, long int *si, long int *sj);

/* End of prototypes */


/*--------------------------------------------------------------------
 * Global Variables
 */

//Strings over the Alphabet Sigma
//char a[] = {'C', 'G', 'T', 'G', 'A', 'A', 'T', 'T', 'C', 'A', 'T'};
//char b[] = {'G', 'A', 'C', 'T', 'T', 'A', 'C'};
char a[m];
char b[n];

//Defines iterators
long int i, j;

/* End of global variables */

/*--------------------------------------------------------------------
 * Function:    main
 */
int main(int argc, char* argv[]) {
    int thread_count = strtol(argv[1], NULL, 10);

    //Generates the values of a
    for (i = 0; i < m; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            a[i] = 'A';
        else if (aux == 2)
            a[i] = 'C';
        else if (aux == 3)
            a[i] = 'G';
        else
            a[i] = 'T';
    }

    //Generates the values of b
    for (i = 0; i < n; i++) {
        int aux = rand() % 4;
        if (aux == 0)
            b[i] = 'A';
        else if (aux == 2)
            b[i] = 'C';
        else if (aux == 3)
            b[i] = 'G';
        else
            b[i] = 'T';
    }

    //Allocates similarity matrix H
    int *H;
    H = malloc((long)m * (long)n * sizeof(long int*));

    //Allocates predecessor matrix P
    int *P;
    P = malloc((long)m * (long)n * sizeof(long int*));

    //Allocates predecessor matrix B
    int *B;
    B = malloc((long)m * (long)n * sizeof(long int*));

    //Start position for backtrack
    long int maxPos = 0;

    //Gets Initial time
    double initialTime = omp_get_wtime();


    //Calculates the similarity Matrix
    long int nDiag = m + n - 1;
    #pragma omp parallel num_threads(thread_count) private(i)
    {
        for (i = 0; i < nDiag; ++i)
        {
            long int nEle = nElement(i);
            #pragma omp for
            for (j = 0; j < nEle; ++j)
            {
                long int si, sj;
                calcIndex(i, j, &si, &sj);
                similarityScore(si, sj, H, P, B, &maxPos);
            }
        }
    }

    backtrack(B, maxPos);

    //Gets final time
    double finalTime = omp_get_wtime();
    printf("\nElapsed time: %f\n", finalTime - initialTime);

    //printf("\nSimilarity Matrix:\n");
    //printMatrix(H);

    //printf("\nPredecessor Matrix:\n");
    //printPredecessorMatrix(P, B);
    //Frees similarity matrixs
    free(H);
    free(P);
    free(B);


    return 0;
}  /* End of main */

/*--------------------------------------------------------------------
 * Function:    nElement
 * Purpose:     Calculate the number of i-diagonal elements
 */
long int nElement(long int i) {
    if (i < m && i < n) {
        //Number of elements in the diagonal is increasing
        return i + 1;
    }
    else if (i < max(m, n)) {
        //Number of elements in the diagonal is stable
        long int min = min(m, n);
        return min;
    }
    else {
        //Number of elements in the diagonal is decreasing
        long int min = min(m, n);
        return 2 * min - i + abs(m - n) - 1;
    }
}

/*--------------------------------------------------------------------
 * Function:    calcElement
 * Purpose:     Calculate the position of (si, sj)-element
 */
void calcIndex(long int i, long int j, long int *si, long int *sj) {
    // Calculate the first element of diagonal
    if (i < n) {
        *si = i;
        *sj = 0;
    } else {
        *si = n - 1;
        *sj = i - n + 1;
    }
    // Calculate the current element of diagonal
    *si -= j;
    *sj += j;
}

/*--------------------------------------------------------------------
 * Function:    SimilarityScore
 * Purpose:     Calculate the maximum Similarity-Score H(i,j)
 */
void similarityScore(long int i, long int j, int* H, int* P, int* B, long int* maxPos) {

    long int up, left, diag;

    //Stores index of element
    long int index = m * i + j;

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
    long int max = NONE;
    long int pred = NONE;
    long int predPos = NONE;
    /* === Matrix ===
     * a[0] ... a[n]
     * b[0]
     * ...
     * b[n]
     *
     * generate 'a' from 'b', if '←' insert e '↑' remove
     * a=GAATTCA
     * b=GACTT-A
     *
     * generate 'b' from 'a', if '←' insert e '↑' remove
     * b=GACTT-A
     * a=GAATTCA
    */
    if (up > max) { //remove letter ↑
        max = up;
        pred = UP;
        predPos = index - m;
    }
    if (left > max) { //insert letter ←
        max = left;
        pred = LEFT;
        predPos = index - 1;
    }
    if (diag > max) { //same letter ↖
        max = diag;
        pred = DIAGONAL;
        predPos = index - m - 1;
    }
    //Inserts the value in the similarity and predecessor matrixes
    H[index] = max;
    P[index] = pred;

    //Check if predecessor position is valid
    //first line has to go left and first col has to go up
    // if ((i != 0 && j != 0) || i == 0 && pred == LEFT || j == 0 && pred == UP)
    B[index] = predPos;
    // else
    // B[index] = 0;

    //Updates maximum score to be used as seed on backtrack
    if (max > H[*maxPos]) {
        #pragma omp critical
        *maxPos = index;
    }

}  /* End of similarityScore */


/*--------------------------------------------------------------------
 * Function:    matchMissmatchScore
 * Purpose:     Similarity function on the alphabet for match/missmatch
 */
long int matchMissmatchScore(long int i, long int j) {
    if (a[j] == b[i]) {
        return matchScore;
    }
    else
        return missmatchScore;
}  /* End of matchMissmatchScore */

/*--------------------------------------------------------------------
 * Function:    backtrack
 * Purpose:     Modify matrix to print, path change from value to PATH
 */
void backtrack(int* B, long int maxPos) {
    //hold maxPos value
    long int predPos;

    //backtrack from maxPos to startPos = 0
    do {
        predPos = B[maxPos];
        B[maxPos] = PATH;
        maxPos = predPos;
    } while (B[maxPos] != 0);

    //sets start '0' to '-1'
    B[maxPos] = PATH;

}  /* End of backtrack */

/*--------------------------------------------------------------------
 * Function:    printMatrix
 * Purpose:     Prlong int Matrix
 */
void printMatrix(int* matrix) {
    long int i, j;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            printf("%d\t", matrix[m * i + j]);
        }
        printf("\n");
    }

}  /* End of printMatrix */

/*--------------------------------------------------------------------
 * Function:    printPredecessorMatrix
 * Purpose:     Prlong int predecessor matrix
 */
void printPredecessorMatrix(int* matrix, int* B) {
    long int i, j, index, path;
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            index = m * i + j;
            path = B[index] + 1;
            if (matrix[index] == UP) {
                if (path) printf("↑ ");
                else printf(BOLDRED "↑ " RESET);
            }
            else if (matrix[index] == LEFT) {
                if (path) printf("← ");
                else printf(BOLDRED "← " RESET);
            }
            else if (matrix[index] == DIAGONAL) {
                if (path) printf("↖ ");
                else printf(BOLDRED "↖ " RESET);
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
    for (i = 0; i < n; i++) { //Lines
        for (j = 0; j < m; j++) {
            printf("%d\t", m * i + j);
        }
        printf("\n");
    }
}  /* End of printPositionMatrix */



/*--------------------------------------------------------------------
 * External References:
 * http://vlab.amrita.edu/?sub=3&brch=274&sim=1433&cnt=1
 * http://pt.slideshare.net/avrilcoghlan/the-smith-waterman-algorithm
 * http://baba.sourceforge.net/
 */