# ff4
Fast f4 algorithm to find groebner basis

Benchmarks divided into internal and external benchmarks. To run internal benchmarks, you need to just clone the repository and build cmake.
If you want to run external benchmarks, and check results by yourself, first you need to download http://nauotit.github.io/openf4/

Current benchmarks:
| Name of test | Time |
| ------------- | ------------- |
 f4_cyclic4  | 7.075 milliseconds. | 
| GroebnerBasisLibF4_cyclic4  | 1652.4 milliseconds. | 
| openf4_cyclic4  | 548.274 milliseconds. | 
| f4_katsura4  | 4.238 milliseconds. | 
| GroebnerBasisLibF4_katsura4  | 2116.03 milliseconds. | 
| openf4_katsura4  | 476.342 milliseconds. | 
| f4_sym3-3  | 126.436 milliseconds. | 
| GroebnerBasisLibF4_sym3-3  | 2077.09 milliseconds. | 
| openf4_sym3-3  | 478.124 milliseconds. | 
| f4_cyclic-5  | 370.958 milliseconds. | 
| openf4_cyclic-5  | 1419.53 milliseconds. | 
| f4_cyclic-6  | 3.039 milliseconds. | 
| openf4_cyclic-6  | 7.071 milliseconds. | 
| f4_cyclic-7  | 73.049 milliseconds. | 
| openf4_cyclic-7  | 231.593 milliseconds. | 
| f4_katsura-9  | 126.555 milliseconds. | 
| openf4_katsura9  | 106.432 milliseconds. | 
| f4_katsura-10  | 652.717 milliseconds. | 
| openf4_katsura10  | 585.24 milliseconds. | 
| f4_katsura-11  | 3514.09 milliseconds. | 
| openf4_katsura11  | 3520.95 milliseconds. | 
| f4_katsura-12  | 21682.2 milliseconds. | 
| openf4_katsura12  | 22284.5 milliseconds. | 
