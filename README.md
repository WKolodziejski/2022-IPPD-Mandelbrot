<p align="center">
  <img src="/img/mandelbrot.png" width="400">
</p>

# Compile:
````
mpic++ main.cpp -o main -fopenmp
````

or

````
make
````

# Run:
````
mpirun --hostfile hosts -np <numero-processos> ./main <tamanho-imagem> <max-iterations>
````