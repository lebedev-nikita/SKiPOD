#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {

    if (argc != 3) {
        printf("Неверное чило параметров\nОжидается: nx, numasks, rank\n");
        return -1;
    }

    int nx = atoi(argv[1]);
    int numtasks = atoi(argv[2]);
    // int rank = atoi(argv[3]);


    printf("start\tlast\tnrows\n");
    for (int rank = 0; rank < numtasks; rank++) {
        int startrow = (rank * nx) / numtasks;          ////
        int lastrow = ( (rank + 1)*nx / numtasks ) - 1; // TODO: здесь может не поделиться 
        int nrows = lastrow - startrow + 1;              ////
        printf("%d\t%d\t%d\n", startrow, lastrow, nrows);
    }
    return 0;
}