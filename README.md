# Smith-Waterman Algorithm in OpenMP
A highly efficient and simple implementation of the Smith-Waterman algorithm for local alignment of sequences.


## Sample Output
![Demo Doccou alpha](https://github.com/TheFighters/Smith-Waterman/blob/master/Media/sampleOutput.png?raw=true)


## How to Use
### Serial Version
* Compile the project
    
    ```bash
    gcc serial_smithW.c -o serial_smithW -fopenmp -DDEBUG
    ```
* Run the program:
    
    ```bash
    ./serial_smithW <number_of_col> <number_of_rows>
    ```

### Parallel Version
* Compile the project
    
    ```bash
    gcc omp_smithW.c -o omp_smithW -fopenmp -DDEBUG
    ```
* Run the program:
    
    ```bash
    ./omp_smithW <number_of_threads> <number_of_col> <number_of_rows>
    ```
    
## Bug Reports & Feature Requests
You can help by reporting bugs, suggesting features, reviewing feature specifications or just by sharing your opinion.
Use [GitHub Issues](https://github.com/TheFighters/Smith-Waterman/issues) for all of that.

## To do
* Display the result in a more user-friendly way.
* Make speedup tests to show how OpenMP accelerated the code.
* Add documentation about the diagonal approach used to implement the parallel algorithm.
    
    
## Contributing
1. Fork the project.
2. Create a branch for your new feature.
3. Test your code.
5. Submit a pull request.

All pull requests are welcome !

## Authors
This project was develloped by [Daniel Holanda](https://github.com/danielholanda/), [Hanoch Griner](https://github.com/eugriner) and [Taynara Pinheiro](https://github.com/taypi).

## License
This project uses the MIT license. See [LICENSE](https://github.com/TheFighters/Smith-Waterman/blob/master/LICENSE) for more details.
