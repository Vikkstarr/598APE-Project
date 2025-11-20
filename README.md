# Gerbil Recreation: A k-mer Counter

This is a recreation of Gerbil, a fast and memory efficient k-mer counter. Counting substrings in genome sequences of length k is required in many bioinformatics problems. We attempt to recreate some of the memory and CPU performance optimizations covered by the original Gerbil paper [here](https://almob.biomedcentral.com/articles/10.1186/s13015-017-0097-9).

## Installation:

Download the source files:

``` git clone https://github.com/Vikkstarr/598APE-Project ```

``` cd src ```


## Compilation:

You can compile the k-mer counting pipeline by running the following command:

``` g++ -std=c++17 -pthread -o pipeline pipeline.cpp Hasher.cpp -O3 ```


## Usage:

You can use the k-mer counting pipeline by running the following command:

Run with an input FASTA file:

``` ./pipeline <input_path> <k> <m> <num_threads> ```

eg. ``` ./pipeline test.fasta 6 5 8 ```



<br>

Run with a randomly-generated file of a specified size:

``` ./pipeline <input_size> <k> <m> <num_threads> ```

eg. ``` ./pipeline 5000000 6 5 8 ```  

---
Citation:

> Marius Erbert, Steffen Rechner, and Matthias Müller-Hannemann, Gerbil: A fast and memory-efficient k-mer counter with GPU-support, Algorithms for Molecular Biology (2017) 12:9, open access.
