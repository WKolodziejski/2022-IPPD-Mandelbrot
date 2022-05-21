#include <iostream>
#include <complex>
#include <string>
#include <mpi.h>

using std::string;
using std::cout;

constexpr auto image_size = 1000, max_iteration = 100;
float offset_x = 1.0f;
float offset_y = 0.5f;
int world_size;

typedef enum {
    DATA_TAG, TERM_TAG, RESULT_TAG
} Tags;

inline int offset2(int x, int y) {
    return x * 3 + y;
}

void worker() {
    MPI_Status status;
    int row;

    int *output = (int *) malloc(image_size * 3 * sizeof(int));

    MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    while (status.MPI_TAG == DATA_TAG) {
        #pragma omp parallel default(none) shared(output, status, row, offset_x, offset_y)
        {
            #pragma omp for
            for (auto column = 0; column < image_size; ++column) {
                std::complex<float> z, c = {
                        (float) column * offset_x / image_size - offset_y,
                        (float) row * offset_x / image_size - (offset_y + 0.5f)
                };

                int i = 0;
                while (abs(z) < 2 && ++i < max_iteration)
                    z = pow(z, 2) + c;

                double t = (double) i / (double) max_iteration;

                int r = (int) (9 * (1 - t) * t * t * t * 255);
                int g = (int) (15 * (1 - t) * (1 - t) * t * t * 255);
                int b = (int) (8.5 * (1 - t) * (1 - t) * (1 - t) * t * 255);

                output[offset2(column, 0)] = r;
                output[offset2(column, 1)] = g;
                output[offset2(column, 2)] = b;
            }
        }
        MPI_Send(&output[0], image_size * 3, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    free(output);
}

void leader() {
    MPI_Status status;
    int count = 0;
    int row = 0;
    int rank;

    for (int i = 1; i < world_size; i++) {
        MPI_Send(&row, 1, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD);
        count++;
        row++;
    }

    int *output = (int *) malloc(image_size * 3 * sizeof(int));
    FILE *file = fopen("mandelbrot.ppm", "w");

    fprintf(file, "P3\n");
    fprintf(file, "%d %d\n", image_size, image_size);
    fprintf(file, "255\n");

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    do {
        MPI_Recv(&output[0], image_size * 3, MPI_INT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
        count--;

        if (row < image_size) {
            MPI_Send(&row, 1, MPI_INT, status.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD);
            row++;
            count++;
        } else {
            MPI_Send(&row, 1, MPI_INT, status.MPI_SOURCE, TERM_TAG, MPI_COMM_WORLD);
        }

        for (auto column = 0; column < image_size; ++column) {
            fprintf(file, "%d\t", output[offset2(column, 0)]);
            fprintf(file, "%d\t", output[offset2(column, 1)]);
            fprintf(file, "%d\t", output[offset2(column, 2)]);
            fprintf(file, "\t");
        }

        fprintf(file, "\n");

    } while (count > 0);

    free(output);
    fclose(file);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    printf("Processor %s, rank %d out of %d processors\n",
           processor_name, world_rank, world_size);

    if (world_rank == 0) {
        leader();
    } else {
        worker();
    }

    MPI_Finalize();
}

//compile: mpic++ main.cpp -o main
//run: mpirun --hostfile hosts -np 0 ./main
