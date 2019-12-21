#include "mpi_fdtd-2d.h"
double bench_t_start, bench_t_end;
int numtasks, rank;
int startrow, lastrow, nrows;

static
double rtclock()
{
        struct timeval Tp;
        int stat;
        stat = gettimeofday (&Tp, NULL);
        if (stat != 0)
            printf ("Error return from gettimeofday: %d", stat);
        return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
    bench_t_start = rtclock ();
}

void bench_timer_stop()
{
    bench_t_end = rtclock ();
}

void bench_timer_print()
{
    printf ("Time in seconds = %0.6lf\n", bench_t_end - bench_t_start);
}

/* Каждый процесс работает только со своей частью матрицы. */
static
void init_array (int tmax, int nx, int ny,
                 float ex[nrows][ny],   // ex не нужны теневые грани
                 float ey[nrows+2][ny], // ey и hz нужны теневые грани
                 float hz[nrows+2][ny], // 
                 float _fict_[tmax]
    )
{
    int i, j;
    for (i = 0; i < tmax; i++)
        _fict_[i] = (float) i;
    for (i = 1; i <= nrows; i++)
        for (j = 0; j < ny; j++)
        {
            ex[i-1][j] = ( (float) (startrow + i - 1) * (j + 1) ) / nx;
            ey[i][j]   = ( (float) (startrow + i - 1) * (j + 2) ) / ny;
            hz[i][j]   = ( (float) (startrow + i - 1) * (j + 3) ) / nx;
        }
    for (j = 0; j < ny; j++) // инициализируем теневые грани
    {
        ey[0][j]         = ( (float) (startrow + 0 - 1) * (j + 2) ) / ny;
        hz[0][j]         = ( (float) (startrow + 0 - 1) * (j + 3) ) / nx;
        ey[nrows+1][j]   = ( (float) (startrow + nrows+1 - 1) * (j + 2) ) / ny;
        hz[nrows+1][j]   = ( (float) (startrow + nrows+1 - 1) * (j + 3) ) / nx;
    }
}

// основные вычисления происходят здесь
static
void kernel_fdtd_2d (int tmax, int nx, int ny,
                     float ex[nrows][ny], float ey[nrows+2][ny], float hz[nrows+2][ny], float _fict_[tmax]
                    )
{
    MPI_Request req[4];
    MPI_Status status[4];
    int t, i, j;
    int nExchanges, n;
    int II, shift;

    for(t = 0; t < tmax; t++)
    {

        // начат участок ey
        if (rank == 0)
            for (j = 0; j < ny; j++)
                ey[1][j] = _fict_[t];
        else
            for (j = 0; j < ny; j++)
                ey[1][j] = ey[1][j] - 0.5f*(hz[1][j]-hz[0][j]);

        for (i = 2; i <= nrows; i++) {
            for (j = 0; j < ny; j++)    
                ey[i][j] = ey[i][j] - 0.5f*(hz[i][j]-hz[i-1][j]); 
        }
        // завершен участок ey
        // начат участок ex
        for (i = 0; i < nrows; i++) { // ex не имеет теневых граней
            for (j = 1; j < ny; j++)
                ex[i][j] = ex[i][j] - 0.5f*(hz[i+1][j]-hz[i+1][j-1]);
        }
        // завершен участок ex

        // начинаем синхронизировать ey
        if (rank != 0) {
            MPI_Irecv(&ey[0][0],       ny, MPI_FLOAT, rank - 1, 'e'+'y'+1, MPI_COMM_WORLD, &req[0]);
        }

        if (rank != 0) {
            MPI_Isend(&ey[1][0],       ny, MPI_FLOAT, rank - 1, 'e'+'y'+2, MPI_COMM_WORLD, &req[1]);
        }

        if (rank != numtasks - 1) {
            MPI_Isend(&ey[nrows][0],   ny, MPI_FLOAT, rank + 1, 'e'+'y'+1, MPI_COMM_WORLD, &req[2]);
        }

        if (rank != numtasks - 1) {
            MPI_Irecv(&ey[nrows+1][0], ny, MPI_FLOAT, rank + 1, 'e'+'y'+2, MPI_COMM_WORLD, &req[3]);
        }

        II = 4;  
        shift = 0;
        if (numtasks-1 == 0) {II = 0; shift=0;} // если существует всего 1 MPI процесс
        else if (rank == 0) {II = 2; shift = 2;}
        else if (rank == numtasks-1) {II = 2; shift = 0;}
        MPI_Waitall(II, &req[shift], &status[0]);
        // закончили синхронизацию ey
        // начат участок hz
        for (i = 1; i <= nrows - 1; i++) {
            for (j = 0; j < ny - 1; j++)
                hz[i][j] = hz[i][j] - 0.7f*(ex[i-1][j+1] - ex[i-1][j] + ey[i+1][j] - ey[i][j]);
        }

        if (rank != numtasks-1) {
            for (j = 0; j < ny - 1; j++)
                hz[nrows][j] = hz[nrows][j] - 0.7f*(ex[nrows-1][j+1] - ex[nrows-1][j] + ey[nrows+1][j] - ey[i][j]);
        }
        // завершен участок hz

        // начинаем синхронизацию hz
        if (rank != 0) {
            MPI_Irecv(&hz[0][0],       ny, MPI_FLOAT, rank - 1, 'h'+'z'+1, MPI_COMM_WORLD, &req[0]);
        }

        if (rank != 0) {
            MPI_Isend(&hz[1][0],       ny, MPI_FLOAT, rank - 1, 'h'+'z'+2, MPI_COMM_WORLD, &req[1]);
        }

        if (rank != numtasks - 1) {
            MPI_Isend(&hz[nrows][0],   ny, MPI_FLOAT, rank + 1, 'h'+'z'+1, MPI_COMM_WORLD, &req[2]);
        }

        if (rank != numtasks - 1) {
            MPI_Irecv(&hz[nrows+1][0], ny, MPI_FLOAT, rank + 1, 'h'+'z'+2, MPI_COMM_WORLD, &req[3]);
        }

        II = 2;
        shift = 0;
        if (numtasks-1 == 0) {II = 0; shift=0;}  // если существует всего 1 MPI процесс
        else if (rank == 0) {II = 2; shift = 2;}
        else if (rank == numtasks-1) {II = 2; shift = 0;}
        MPI_Waitall(II, &req[shift], &status[0]);
        // завершили синхронизацию hz
    }
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
    int leader_rank = 0;

    int tmax = TMAX;
    int nx = NX;
    int ny = NY;

    startrow = (rank * nx) / numtasks;
    lastrow = ( (rank + 1)*nx / numtasks ) - 1;
    nrows = lastrow - startrow + 1;


    float (*ex)[nrows][ny];     ex = ( float(*)[nrows+2][ny] ) malloc ( (nrows) * ny * sizeof(float));
    float (*ey)[nrows+2][ny];   ey = ( float(*)[nrows+2][ny] ) malloc ( (nrows+2) * ny * sizeof(float));
    float (*hz)[nrows+2][ny];   hz = ( float(*)[nrows][ny]   ) malloc ( (nrows+2) * ny * sizeof(float));
    float (*_fict_)[tmax];   _fict_ = ( float(*)[tmax]       ) malloc (tmax * sizeof(float));
    init_array (tmax, nx, ny, *ex, *ey, *hz, *_fict_);


    if (rank == leader_rank)
        bench_timer_start();

    kernel_fdtd_2d (tmax, nx, ny, *ex, *ey, *hz, *_fict_); // вычисления

    MPI_Barrier(MPI_COMM_WORLD); // ждем, чтобы измерения времени были максимально точными

    if (rank == leader_rank) 
    {
        bench_timer_stop();
        bench_timer_print();
    }

    free( (void*)ex     );
    free( (void*)ey     );
    free( (void*)hz     );
    free( (void*)_fict_ );

    MPI_Finalize();
    return 0;
}
