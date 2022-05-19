#include <iostream>
#include <complex>
#include <string>
#include <omp.h>

using std::string;
using std::cout;

constexpr auto max_row = 30, max_column = 78, max_iteration = 255;

int main() {
    u_int8_t output[max_row][max_column];

    #pragma omp parallel default(none) shared(output, cout, offset_x, offset_y)
    {
        #pragma omp for
        for (auto row = 0; row < max_row; ++row) {

            #pragma omp parallel default(none) shared(output, cout, row, offset_x, offset_y)
            {
                #pragma omp for
                for (auto column = 0; column < max_column; ++column) {
                    std::complex<float> z, c = {
                            (float) column * offset_x / max_column - offset_y,
                            (float) row * offset_x / max_row - (offset_y + 0.5f)
                    };

                    int i = 0;
                    while (abs(z) < 2 && ++i < max_iteration)
                        z = pow(z, 2) + c;

                    double t = (double)i/(double)max_iteration;

                    int r = (int)(9*(1-t)*t*t*t*255);
                    int g = (int)(15*(1-t)*(1-t)*t*t*255);
                    int b =  (int)(8.5*(1-t)*(1-t)*(1-t)*t*255);

                    output[row * max_column + column] = r;
                    output[row][column][1] = g;
                    output[row][column][2] = b;
                }
            }
        }
    }


    /*for (auto & row : output) {
        for (char column : row) {
            std::cout << column;
        }
        std::cout << '\n';
    }*/
    FILE *file = fopen("test.ppm", "w");
    fprintf(file, "P3\n");
    fprintf(file, "%d %d\n", max_row, max_column);
    fprintf(file, "255\n");

    for (auto &row: output) {
        for (int *column: row) {
            //std::cout << column;
            fprintf(file, "%d\t", column[0]);
            fprintf(file, "%d\t", column[1]);
            fprintf(file, "%d\t", column[2]);
            fprintf(file, "\t");
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0;
}
