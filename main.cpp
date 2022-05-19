#include <iostream>
#include <complex>
#include <string>
#include <omp.h>

using std::string;
using std::cout;

constexpr auto max_row = 30, max_column = 78, max_iteration = 255;

int main() {
    u_int8_t output[max_row][max_column];

    #pragma omp parallel default(none) shared(output, cout)
    {
        #pragma omp for
        for (auto row = 0; row < max_row; ++row) {

            #pragma omp parallel default(none) shared(output, cout, row)
            {
                #pragma omp for
                for (auto column = 0; column < max_column; ++column) {
                    std::complex<float> z, c = {
                            (float) column * 2.0f / max_column - 1.5f,
                            (float) row * 2.0f / max_row - 1
                    };

                    int i = 0;
                    while (abs(z) < 2 && ++i < max_iteration)
                        z = pow(z, 2) + c;

                    //output[row][column] = (i == max_iteration ? '#' : '.');
                    //output[row][column] = (i == max_iteration ? '#' : '.');
                    output[row][column] = i;
                    cout << i << ' ';
                    //cout << output[row][column];
                }
            }
            cout << '\n';
        }
    }


    /*for (auto & row : output) {
        for (char column : row) {
            std::cout << column;
        }
        std::cout << '\n';
    }*/
    FILE *file = fopen("test.ppm", "w");

    char* header = "P3 \n 30 78 \n 255 \n";
    cout << std::endl << sizeof(header) << std::endl;
    fwrite(header, sizeof(*header), 1, file);

    for (auto & row : output) {
        for (char column : row) {
            //std::cout << column;
            fwrite(&column, sizeof(column), 1, file);
            fwrite(&column, sizeof(column), 1, file);
            fwrite(&column, sizeof(column), 1, file);
        }
        //std::cout << '\n';
    }
    //fwrite(<>, sizeof(*<>), 1, file);
    fclose(file);

    return 0;
}
