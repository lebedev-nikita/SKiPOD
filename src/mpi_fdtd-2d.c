#include "mpi_fdtd-2d.h"
double bench_t_start, bench_t_end;
int numtasks, rank;                             ////
int startrow, lastrow, nrows;                   ////

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

static
void init_array (int tmax, int nx, int ny,
                 float ex[nrows][ny],
                 float ey[nrows+2][ny],
                 float hz[nrows+2][ny],
                 float _fict_[tmax]
    )                                           ////
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

static
void print_array(int nx,
                 int ny,
                 float ex[nx][ny],
                 float ey[nx][ny],
                 float hz[nx][ny]
                )
{
    int i, j;
    printf("==BEGIN DUMP_ARRAYS==\n");

    printf("begin dump: %s\n", "ex");
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            printf("%5.2f ", ex[i][j]);
        }
        printf("\n");
    }
    printf("end   dump: %s\n", "ex");

    printf("begin dump: %s\n", "ey");
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            printf("%5.2f ", ey[i][j]);
        }
        printf("\n");
    }
    printf("end   dump: %s\n", "ey");

    printf("begin dump: %s\n", "hz");
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            printf("%5.2f ", hz[i][j]);
        }
        printf("\n");
    }
    printf("end   dump: %s\n", "hz");
    
    printf("==END   DUMP_ARRAYS==\n");
}

// основные вычисления
static
void kernel_fdtd_2d (int tmax, int nx, int ny,
                     float ex[nrows][ny], float ey[nrows+2][ny], float hz[nrows+2][ny], float _fict_[tmax]
                    )
{
    int t, i, j;
    for(t = 0; t < tmax; t++)
    {
        if (rank == 0)                                              ////
            for (j = 0; j < ny; j++)
                ey[1][j] = _fict_[t];
        else
            for (j = 0; j < ny; j++)
                ey[1][j] = ey[i][j] - 0.5f*(hz[1][j]-hz[0][j]); ////

        for (i = 2; i <= nrows; i++) {
            for (j = 0; j < ny; j++)    
                ey[i][j] = ey[i][j] - 0.5f*(hz[i][j]-hz[i-1][j]); ////
        }
        for (i = 0; i < nrows; i++) {                              ////
            for (j = 1; j < ny; j++)
                ex[i][j] = ex[i][j] - 0.5f*(hz[i+1][j]-hz[i+1][j-1]);
        }

        for (i = 1; i <= nrows - 1; i++) {                           ////
            for (j = 0; j < ny - 1; j++)
                hz[i][j] = hz[i][j] - 0.7f*(ex[i-1][j+1] - ex[i-1][j] + ey[i+1][j] - ey[i][j]);
        }
        // V
    }
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);                     ////
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);   ////
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);       ////
    MPI_Barrier(MPI_COMM_WORLD);                ////
    int leader_rank = 0;             ////

    int tmax = TMAX;
    int nx = NX;
    int ny = NY;

    startrow = (rank * nx) / numtasks;          ////
    lastrow = ( (rank + 1)*nx / numtasks ) - 1; ////
    nrows = lastrow - startrow + 1;              ////

    // #pragma omp parallel
    //     printf("%d\t %d\n", rank, omp_get_thread_num());
    float (*ex)[nrows][ny];     ex = ( float(*)[nrows+2][ny] ) malloc ( (nrows) * ny * sizeof(float));
    float (*ey)[nrows+2][ny];   ey = ( float(*)[nrows+2][ny] ) malloc ( (nrows+2) * ny * sizeof(float));
    float (*hz)[nrows+2][ny];   hz = ( float(*)[nrows][ny]   ) malloc ( (nrows+2) * ny * sizeof(float));
    float (*_fict_)[tmax];   _fict_ = ( float(*)[tmax]       ) malloc (tmax * sizeof(float));
    init_array (tmax, nx, ny, *ex, *ey, *hz, *_fict_);

    if (rank == leader_rank)
        bench_timer_start();

    kernel_fdtd_2d (tmax, nx, ny, *ex, *ey, *hz, *_fict_); // вычисления

    MPI_Barrier(MPI_COMM_WORLD);                ////
    if (rank == leader_rank) 
    {
        bench_timer_stop();
        bench_timer_print();
        // print_array(nx, ny, *ex, *ey, *hz);
    }

    free( (void*)ex     );
    free( (void*)ey     );
    free( (void*)hz     );
    free( (void*)_fict_ );

    MPI_Finalize();                             ////
    return 0;
}
