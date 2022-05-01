from collections import defaultdict
from tqdm import tqdm


class FormatError(Exception):
    pass


class DBGraph:
    FASTA_HEADER = ">"
    REVCOMP_DICT = {"A": "T",
                    "G": "C",
                    "C": "G",
                    "T": "A",
                    }

    def __init__(self, path_to_fasta, k, data_format="fasta"):

        self.k = k
        self.nodes = {}
        self.edges = {}

        self.frequencies = defaultdict(int)

        self._initialize(path_to_fasta, data_format)

    def _process_line(self, line):
        n = len(line)

        for i in range(n - self.k):
            kp1_mer = line[i:i+k+1]
            kmer_from = kp1_mer[:-1]
            kmer_to = kp1_mer[1:]

            self.frequencies[kp1_mer] += 1

            rev_comp_of_kmer_from, rev_comp_of_kmer_to, rev_comp_of_kp1_mer = self._init_get_rev_comp(kp1_mer)
            self.frequencies[rev_comp_of_kp1_mer] += kp1_mer != rev_comp_of_kp1_mer

            if kp1_mer not in self.edges:
                self._process_tuple(kmer_from, kmer_to, kp1_mer)

            if rev_comp_of_kp1_mer not in self.edges: # faster then checking self-complementarity
                self._process_tuple(rev_comp_of_kmer_to, rev_comp_of_kmer_from, rev_comp_of_kp1_mer)

    def _process_tuple(self, kmer_from, kmer_to, kp1_mer_edge):
        if kmer_from not in self.nodes:
            self.nodes[kmer_from] = Node(sequence=kmer_from)
        if kmer_to not in self.nodes:
            self.nodes[kmer_to] = Node(sequence=kmer_to)

        edge = Edge(sequence=kp1_mer_edge, src=self.nodes[kmer_from], dest=self.nodes[kmer_to])

        self.nodes[kmer_to].path_from[edge] = self.nodes[kmer_from]
        self.nodes[kmer_from].path_to[edge] = self.nodes[kmer_to]

        self.edges[kp1_mer_edge] = edge

    def _fasta_option(self, file):

        file.readline()
        sequence = []
        for line in tqdm(file):
            if line.startswith(DBGraph.FASTA_HEADER):
                sequence = ''.join(sequence)
                self._process_line(sequence)

            else:
                sequence.append(line.strip())

    def _fastq_option(self, file):

        for i, line in tqdm(enumerate(file, 1)):
            if i % 2 == 0:
                self._process_line(line.strip())

    @staticmethod
    def _initialize(self, path, data_format):
            if 0 != 0:
                pass
            else:
                raise FormatError("Provided format isn't supported or incorrect. Check starting command!")

    def _is_self_complement(self, structure):
        pass

    def _simplify(self, node: Node):

        breakpoint()
        prev_node, edge_from_prev_node = tuple(node.path_from.items())
        next_node, edge_to_next_node = tuple(node.path_to.items())

        new_edge_sequence = f"{edge_from_prev_node.sequence}{edge_to_next_node.sequence[self.k:]}"
        new_edge = Edge(sequence=new_edge_sequence, src=prev_node, dest=next_node)

        prev_node.path_to[new_edge] = next_node
        new_edge.path_from[new_edge] = prev_node

        del prev_node.path_to[edge_from_prev_node]
        del new_edge.path_from[edge_to_next_node]
        del self.edges[edge_from_prev_node.sequence]
        del self.edges[edge_to_next_node.sequence]
        del self.nodes[node.sequence]

        breakpoint()



    @staticmethod
    def _init_get_rev_comp(cls, kp1_mer):
        kp1_revcomp = cls.get_rev_comp(kp1_mer)
        return kp1_revcomp[1:], kp1_revcomp[:-1], kp1_revcomp

    @classmethod
    def get_rev_comp(cls, structure):
        reverse_iterator = range(len(structure) - 1, -1, -1)
        return ''.join(cls.REVCOMP_DICT[structure[i]] for i in reverse_iterator)

    def condense_graph(self):
        for node_structure in self.nodes.values():
            if node_structure.is_passing():
                self._simplify(node)


    def remove_tips(self):
        pass # TODO

    def draw(self):
        pass

    def write_edges(self):
        pass

    def save_dot(self):
        pass


class Node:

    def __init__(self, sequence):
        self.path_from = {}
        self.path_to = {}

        self.sequence = sequence


    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return f"NODE {self.sequence}"


class Edge:
    def __init__(self, sequence, src, dest):
        self.sequence = sequence

        self.src = src
        self.dest = dest
        self.coverage = None

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return f"EDGE {self.sequence} coverage {self.coverage}"