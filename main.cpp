#include <iostream>
#include <complex>
#include <string>
#include <map>
#include <mpi.h>

using std::cout;
using std::string;

int image_size = 1000, max_iteration = 100;
float p_x1 = -1.0f;
float p_y1 = -1.0f;
float p_x2 = 1.0f;
float p_y2 = 1.0f;
/*float offset_x = 1.0f;
float offset_y = 0.5f;*/
int world_size;

typedef enum {
    DATA_TAG,
    TERM_TAG,
    RESULT_TAG
} Tags;

inline int offset(int x, int y) {
    return x * 3 + y;
}

void worker() {
    MPI_Status status;
    int row;

    int total_elements = (image_size * 3 + 1);
    int *buffer = (int *) malloc(total_elements * sizeof(int));
    int *output = buffer + 1;

    MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

    while (status.MPI_TAG == DATA_TAG) {
        buffer[0] = row;
        #pragma omp parallel default(none) shared(buffer, output, status, row, max_iteration, image_size, p_x1, p_y1, p_x2, p_y2)
        {
            #pragma omp for
            for (auto column = 0; column < image_size; ++column) {
                std::complex<float> z, c = {p_x1 + ((float) column / (float) image_size) * (p_x2 - p_x1),
                                            p_y1 + ((float) row / (float) image_size) * (p_y2 - p_y1)};
                /*std::complex<float> z, c = {
                                           (float)column * offset_x / image_size - offset_y,
                                           (float)row * offset_x / image_size - (offset_y + 0.5f)};*/

                /*   X1 + (x / WIDTH) * (X2 - X1),
                     Y1 + (y / HEIGHT) * (Y2 - Y1)*/

                int i = 0;
                while (abs(z) < 2 && ++i < max_iteration)
                    z = pow(z, 2) + c;

                double t = (double) i / (double) max_iteration;

                int r = (int) (9 * (1 - t) * t * t * t * 255);
                int g = (int) (15 * (1 - t) * (1 - t) * t * t * 255);
                int b = (int) (8.5 * (1 - t) * (1 - t) * (1 - t) * t * 255);

                output[offset(column, 0)] = r;
                output[offset(column, 1)] = g;
                output[offset(column, 2)] = b;
            }
        }
        MPI_Send(&buffer[0], total_elements, MPI_INT, 0, RESULT_TAG, MPI_COMM_WORLD);
        MPI_Recv(&row, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    }

    free(buffer);
}

void printRow(FILE *file, int *output) {
    for (auto column = 0; column < image_size; ++column) {
        fprintf(file, "%d\t", output[offset(column, 0)]);
        fprintf(file, "%d\t", output[offset(column, 1)]);
        fprintf(file, "%d\t", output[offset(column, 2)]);
        fprintf(file, "\t");
    }
}

void leader() {
    MPI_Status status;
    int count = 0;
    int row = 0;
    //int rank;
    std::map<int, int *> row_buffer;
    int received_row;
    int next_row = 0;

    for (int i = 1; i < world_size; i++) {
        MPI_Send(&row, 1, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD);
        count++;
        row++;
    }

    int total_elements = (image_size * 3 + 1);
    int *buffer = (int *) malloc(total_elements * sizeof(int));
    int *output = buffer + 1;

    FILE *file = fopen("mandelbrot.ppm", "w");

    fprintf(file, "P3\n");
    fprintf(file, "%d %d\n", image_size, image_size);
    fprintf(file, "255\n");

    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    do {
        MPI_Recv(&buffer[0], total_elements, MPI_INT, MPI_ANY_SOURCE, RESULT_TAG, MPI_COMM_WORLD, &status);
        received_row = buffer[0];
        count--;

        if (row < image_size) {
            MPI_Send(&row, 1, MPI_INT, status.MPI_SOURCE, DATA_TAG, MPI_COMM_WORLD);
            row++;
            count++;
        } else {
            MPI_Send(&row, 1, MPI_INT, status.MPI_SOURCE, TERM_TAG, MPI_COMM_WORLD);
        }

        if (received_row != next_row) {
            row_buffer.insert({received_row, buffer});
            buffer = (int *) malloc(total_elements * sizeof(int));
            output = buffer + 1;
        } else {
            printRow(file, output);
            ++next_row;

            for (auto it = row_buffer.begin(); it != row_buffer.end();/**/) {
                if (it->first == next_row) {
                    printRow(file, it->second + 1);
                    free(it->second);
                    row_buffer.erase(it++);

                    ++next_row;
                } else {
                    break;
                }
            }
        }

        fprintf(file, "\n");

    } while (count > 0);

    free(buffer);
    fclose(file);
}

int main(int argc, char **argv) {
    image_size = (argc < 2) ? image_size : atoi(argv[1]);
    max_iteration = (argc < 3) ? max_iteration : atoi(argv[2]);
    p_x1 = (argc < 4) ? p_x1 : atof(argv[3]);
    p_y1 = (argc < 5) ? p_y1 : atof(argv[4]);
    p_x2 = (argc < 6) ? p_x2 : atof(argv[5]);
    p_y2 = (argc < 7) ? p_y2 : atof(argv[6]);
    /*offset_x        = (argc < 4) ? 1.0f : atof(argv[3]);
    offset_y        = (argc < 5) ? 0.5f : atof(argv[4]);*/

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

// compile: mpic++ main.cpp -o main -fopenmp
// run: mpirun --hostfile hosts -np <numero-processos> ./main <tamanho-imagem> <max-iterations>
