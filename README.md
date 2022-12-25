# COMP-5704: Parallel Algorithms and Applications in Data Science
This repository is the final project for course <a href="https://www.dehne.ca/comp5704">COMP-5704</a>, supervised by Prof F.Dehne.

## Introduction
String matching is a fundamental problem in computer science, with applications in fields such as bioinformatics, search engines, and large database searches. While sequence algorithms are effective at solving this problem, their performance can be limited by the linear increase in the size of the matching string, leading to long search times for large datasets. In recent years, there has been a trend toward using modern CPUs with multiple cores and GPUs for data processing, and parallelization has become an important consideration in algorithm design. Given the data structure of string matching, it can potentially be highly amenable to parallel programming, allowing for efficient solutions to large-scale matching problems on large-scale systems.

The goal of this project was to investigate the parallel implementation of traditional sequence string matching algorithms, implementing their parallel versions using Cilk and Cuda, and comparing the benefits of different string matching algorithms in different parallel architectures. I also examined the potential applications and limitations of these algorithms on a distributed architecture. In this project, I implemented seven different string-matching algorithms in CPU parallel and five in GPU. For algorithms that are difficult to parallelize algorithmically, I partitioned the index to perform parallel operations. Some of the algorithms also used bit parallel or SSE instructions to improve performance. I analyzed the speedup obtained by different algorithms and the worst case speedup obtained.

## Usage
### CPU Version:
Install Cilk, follow instruction from<a href="https://www.opencilk.org/doc/users-guide/install/">opencilk.org.</a>
Compile with SSE 4.2. for example:

```sh
/opt/opencilk/bin/clang++ -fopencilk -O3 StringMatch.cpp -o StringMatch -msse4.2
```
### GPU Version:
Install Cuda

Compile with nvcc:

```sh
nvcc ParallelStringMatch.cu
```


## License

[MIT](LICENSE) ©

