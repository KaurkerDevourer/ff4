# ff4
Fast f4 algorithm to find groebner basis

If you want to build current benchmarks, and check results by yourself, first you need to download http://nauotit.github.io/openf4/

Current benchmarks:
| Name of test | Time |
| ------------- | ------------- |
| f4_cyclic4  | 7.144 milliseconds. | 
| GroebnerBasisLibF4_cyclic4  | 946.379 milliseconds. | 
| openf4_cyclic4  | 351.262 milliseconds. | 
| f4_katsura4  | 60.299 milliseconds. | 
| GroebnerBasisLibF4_katsura4  | 1348.54 milliseconds. | 
| openf4_katsura4  | 311.183 milliseconds. | 
| f4_sym3-3  | 56.29 milliseconds. | 
| GroebnerBasisLibF4_sym3-3  | 1253.38 milliseconds. | 
| openf4_sym3-3  | 303.512 milliseconds. | 
| f4_cyclic-5  | 211.601 milliseconds. | 
| openf4_cyclic-5  | 1007.45 milliseconds. | 
| f4_cyclic-6  | 2.088 milliseconds. | 
| openf4_cyclic-6  | 5.333 milliseconds. | 
| f4_cyclic-7  | 59.085 milliseconds. | 
| openf4_cyclic-7  | 191.26 milliseconds. | 
| f4_katsura-9  | 100.849 milliseconds. | 
| openf4_katsura9  | 91.531 milliseconds. | 
| f4_katsura-10  | 539.973 milliseconds. | 
| openf4_katsura10  | 530.714 milliseconds. | 
| f4_katsura-11  | 3158.75 milliseconds. | 
| openf4_katsura11  | 3238.93 milliseconds. | 
| f4_katsura-12  | 19624.4 milliseconds. | 
| openf4_katsura12  | 21597.2 milliseconds. | 

Make sure, to run the binary with the highest priority. F.e. sudo nice -n -20 ./internal_benchmark
