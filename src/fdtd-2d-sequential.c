#include "fdtd-2d.h"
double bench_t_start, bench_t_end;

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
                 float ex[nx][ny],
                 float ey[nx][ny],
                 float hz[nx][ny],
                 float _fict_[tmax]
    )
{
    int i, j;
    for (i = 0; i < tmax; i++)
        _fict_[i] = (float) i;
    for (i = 0; i < nx; i++)
        for (j = 0; j < ny; j++)
        {
            ex[i][j] = ((float) i*(j+1)) / nx;
            ey[i][j] = ((float) i*(j+2)) / ny;
            hz[i][j] = ((float) i*(j+3)) / nx;
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
                     float ex[nx][ny], float ey[nx][ny], float hz[nx][ny], float _fict_[tmax]
                    )
{
    int t, i, j;
    for(t = 0; t < tmax; t++)
    {
        for (j = 0; j < ny; j++)
            ey[0][j] = _fict_[t];

        for (i = 1; i < nx; i++)
            for (j = 0; j < ny; j++)    
                ey[i][j] = ey[i][j] - 0.5f*(hz[i][j]-hz[i-1][j]);
        //1
        for (i = 0; i < nx; i++)
            for (j = 1; j < ny; j++)
                ex[i][j] = ex[i][j] - 0.5f*(hz[i][j]-hz[i][j-1]);
        //2
        for (i = 0; i < nx - 1; i++)
            for (j = 0; j < ny - 1; j++)
                hz[i][j] = hz[i][j] - 0.7f*(ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
        //3
    }
}


int main(int argc, char** argv)
{
    int tmax = TMAX;
    int nx = NX;
    int ny = NY;

    float (*ex)[nx][ny];       ex = ( float(*)[nx][ny] ) malloc (nx * ny * sizeof(float));
    float (*ey)[nx][ny];       ey = ( float(*)[nx][ny] ) malloc (nx * ny * sizeof(float));
    float (*hz)[nx][ny];       hz = ( float(*)[nx][ny] ) malloc (nx * ny * sizeof(float));
    float (*_fict_)[tmax]; _fict_ = ( float(*)[tmax]   ) malloc (tmax * sizeof(float));

    init_array (tmax, nx, ny, *ex, *ey, *hz, *_fict_);

    bench_timer_start();

    kernel_fdtd_2d (tmax, nx, ny, *ex, *ey, *hz, *_fict_);

    bench_timer_stop();
    bench_timer_print();
    print_array(nx, ny, *ex, *ey, *hz);

    free( (void*)ex     );
    free( (void*)ey     );
    free( (void*)hz     );
    free( (void*)_fict_ );

    return 0;
}
