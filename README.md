# ff4
Fast f4 algorithm to find groebner basis

If you want to build current benchmarks, and check results by yourself, first you need to download http://nauotit.github.io/openf4/

Current benchmarks:
| Name of test | Time |
| ------------- | ------------- |
| f4_cyclic4  | 8.101 milliseconds. | 
| GroebnerBasisLibF4_cyclic4  | 944.359 milliseconds. | 
| openf4_cyclic4  | 353.999 milliseconds. | 
| f4_katsura4  | 58.842 milliseconds. | 
| GroebnerBasisLibF4_katsura4  | 1370.72 milliseconds. | 
| openf4_katsura4  | 317.737 milliseconds. | 
| f4_sym3-3  | 58.884 milliseconds. | 
| GroebnerBasisLibF4_sym3-3  | 1269.18 milliseconds. | 
| openf4_sym3-3  | 312.628 milliseconds. | 
| f4_cyclic-5  | 207.443 milliseconds. | 
| openf4_cyclic-5  | 1028 milliseconds. | 
| f4_cyclic-6  | 2.056 milliseconds. | 
| openf4_cyclic-6  | 5.254 milliseconds. | 
| f4_cyclic-7  | 58.335 milliseconds. | 
| openf4_cyclic-7  | 189.112 milliseconds. | 
| f4_katsura-9  | 101.037 milliseconds. | 
| openf4_katsura9  | 91.491 milliseconds. | 
| f4_katsura-10  | 540.266 milliseconds. | 
| openf4_katsura10  | 522.516 milliseconds. | 
| f4_katsura-11  | 3187.47 milliseconds. | 
| openf4_katsura11  | 3403.97 milliseconds. | 
| f4_katsura-12  | 19654.2 milliseconds. | 
| openf4_katsura12  | 21175.7 milliseconds. | 

Make sure, to run the binary with the highest priority. F.e. sudo top -n -20 ./internal_benchmark
