# BODBuilder - Basic-Operational De-Bruijn graph Builder

## Motivation and quick explanation
This small tool has been created as a part of NGS data analysis homework in Bioinformatics Institute :)


BODBuilder provides pipeline for constructing De-Bruijn according to its mathematical definition and also utilizes heuristics to simplify it and prune from presumably wrong connections. 

## Typical BODBuilder workflow:
1. Takes a single or multiple files with supported data formats (no need to explicitly specify format).
2. Builds primary oriented multigraph **G** - creates connections in a graph using nodes with size `k` and primary edges with sizes `k+1` between them. Each kmer is processed simultaneously with own reverse-complement one.
3. Simplifies graph until convergence - contracts edges and creates minor graph **G\`** of a graph **G**. Its property - no passing vertices <img src="https://render.githubusercontent.com/render/math?math=(deg(v)_{out} = deg(v)_{in} = 1)">.
4. Removes tips from graph (low-covered edges - user-specified threshold).
5. Tries to remove bulges from graph by removing every badly covered edge.
6. Stores graph in `.dot` format.
7. If option `--draw` is provided, draws graph in `.png` format.


Final graph is guaranteed to have only high-covered edges related to mean edge coverage investigated after condensing and removing tips and low-covered edges. Moreover, final graph doesn\'t have passing vertices.

## Usage
1. Clone repository: 
```{bash}
git clone https://github.com/Abusagit/BODBuilder.git
```
2. Install requirements (_if needed_):
```{bash}
pip install numpy tqdm networkx
```

3. Typical command: 

Keep `--draw` to obtain picture in `.png` format (_I assume you added tool directory to PATH or provided absolute path_):
```{bash}
py <path_to_repository>/build_graph.py -i <input_sequence_file> [<another_sequence_file>] \
-o <output_directory> \
-k <kmer sise - odd> \
-b <lower bound of coverage for edges removing> \
--draw 
```

**Example of graph with `k=55`**
![graph_planar_view_k_55](https://user-images.githubusercontent.com/67659154/168481757-0d69cd66-41b2-442c-8d69-7817b79a623f.png)



## Arguments
| Argument | Description |
| ----------- | ----------- |
|`-h, --help`|Show this help message and exit |
|`-i, --input [fasta\|fastq\|fna\|fa]` |Path to your file(s) with data for graph |
|`-o, --outdir [./]` | Output directory |
|`-k , --kmer-size` | Size of a kmer to be used for graph building | 
|`-b , --bad_cov [100]` | Threshold coverage for clipping edge |
|`--draw` | Stay with this option to get picture of a final graph | 
|`--force` | Force override dir if it is\`nt empty | 
