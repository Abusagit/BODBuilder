from collections import defaultdict
from tqdm import tqdm
import logging
import numpy as np
import networkx as nx
import pygraphviz as pgv
from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

import os

logger = logging.getLogger(__name__)


class FormatError(Exception):
    pass


class Node:

    def __init__(self, sequence, rc_sequence):
        self.path_from = {}
        self.path_to = {}

        self.sequence = sequence
        self.is_self_complement = rc_sequence == sequence

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return f"|NODE {self.sequence} | deg_in {self.deg_in} | deg_out {self.deg_out}|"

    @property
    def deg_in(self):
        return len(self.path_from)

    @property
    def deg_out(self):
        return len(self.path_to)

    def is_passing(self):
        return self.deg_in == self.deg_out == 1

    def is_alone(self):
        return self.deg_in == self.deg_out == 0


class Edge:

    def __init__(self, sequence: str, src: Node, dest: Node, kp1: int, rev_comp_sequence: str, coverage=None):
        self.kp1 = kp1
        self.sequence = sequence
        self.src = src
        self.dest = dest

        self.coverage = coverage

        self.is_self_complement = rev_comp_sequence == sequence

        # self.rev_comp_edge = None

    def __hash__(self):
        return hash(self.sequence)

    # @property
    # def coverage(self):
    #     rev_com_coverage = 0 if not self.rev_comp_edge else np.mean(self.rev_comp_edge._coverage)
    #
    #     return np.mean(self._coverage) + rev_com_coverage

    def __len__(self):
        return len(self.sequence) - self.kp1 + 1

    def __repr__(self):
        return f"EDGE {self.sequence}_coverage_{self.coverage}"


class DBGraph:
    FASTA_HEADER = ">"
    REVCOMP_DICT = {"A": "T",
                    "G": "C",
                    "C": "G",
                    "T": "A",
                    }

    STARS = '*' * 100

    def __init__(self, path_to_fasta, k, outdir: str, ratio, data_format):

        self.k = k
        self.nodes = {}
        self.edges = {}

        self.frequencies = defaultdict(int)
        self.rev_com_dict = {}

        self.outdir = outdir
        self.ratio = ratio

        logger.info(f"{DBGraph.STARS}\n")
        logger.info(f"Initializing graph building: reading file {path_to_fasta}")
        self._initialize(path_to_fasta, data_format)
        logger.info(f"Done!\n{DBGraph.STARS}\n\n{DBGraph.STARS}")

        self.coverages = None
        self.coverage_mean = None
        self.coverage_std = None

    def complete_graph_building(self):
        logger.info("Starting graph simplification...")
        logger.info("Condensing graph...")

        self.condense_graph()

        logger.info("Exploring edge coverage distribution...")
        self.assign_coverage_to_edges()  # ???????????????

        self.coverages = np.array([e.coverage for e in self.edges.values()])
        self.coverage_mean = np.mean(self.coverages)
        self.coverage_std = np.std(self.coverages)

        logger.debug(
            f"Mean coverage is {self.coverage_mean}, std {self.coverage_std}, median is {np.median(self.coverages)}"
        )

        # self.draw_and_save_graph()
        # breakpoint()

        logger.info("Performing tips removing procedure...")
        self.remove_tips(condition=self.bad_edge_coverage)
        #
        # self.draw_and_save_graph()
        # breakpoint()

        logger.info("Removing bulges...")
        self.remove_bulges()

        logger.info("Assigning coverages to remaining edges")
        self.assign_coverage_to_edges()
        logger.info(f"Ended graph cleaning! Ready for saving and drawing...\n{DBGraph.STARS}\n")

        logger.debug(f"{np.mean(np.array([e.coverage for e in self.edges.values()]))}, \
        {np.median(np.array([e.coverage for e in self.edges.values()]))} \
        {np.std(np.array([e.coverage for e in self.edges.values()]))}")

        self.draw_and_save_graph()
        self.write_edges()

    def bad_edge_coverage(self, edge: Edge):
        return edge.coverage <= self.ratio

    def _process_kmers_from_sequence(self, line):
        n = len(line)

        for i in range(n - self.k):
            kp1_mer = line[i:i + self.k + 1]
            kmer_from = kp1_mer[:-1]
            kmer_to = kp1_mer[1:]

            rev_comp_of_kmer_from, rev_comp_of_kmer_to, rev_comp_of_kp1_mer = self._init_get_rev_comp(kp1_mer)

            if kp1_mer not in self.edges:
                self._process_tuple(kmer_from, kmer_to, kp1_mer, rc_of_kp_1=rev_comp_of_kp1_mer,
                                    rc_of_from=rev_comp_of_kmer_from, rc_of_to=rev_comp_of_kmer_to)

            # self.edges[kp1_mer]._coverage[0] += 1

            self.frequencies[kp1_mer] += 1

            if not self.edges[kp1_mer].is_self_complement and rev_comp_of_kp1_mer not in self.edges:

                self._process_tuple(rev_comp_of_kmer_to, rev_comp_of_kmer_from, rev_comp_of_kp1_mer,
                                    rc_of_kp_1=kp1_mer, rc_of_from=kmer_to, rc_of_to=kmer_from)

                # if all((self.edges[kp1_mer].rev_comp_edge is not None,
                #         self.edges[rev_comp_of_kp1_mer].rev_comp_edge is not None)):  # assign reverse complement
                #     # edges only once
                #
                #     self.edges[kp1_mer].rev_comp_edge = self.edges[rev_comp_of_kp1_mer]
                #     self.edges[rev_comp_of_kp1_mer] = self.edges[kp1_mer]

                for forward_comp, reverse_comp in zip(
                        (kmer_from, kmer_to, kp1_mer),
                        (rev_comp_of_kmer_to, rev_comp_of_kmer_from, rev_comp_of_kp1_mer)
                ):
                    self.rev_com_dict[forward_comp] = reverse_comp
                    self.rev_com_dict[reverse_comp] = forward_comp

    def _process_tuple(self, kmer_from: str, kmer_to: str, kp1_mer_edge: str,
                       rc_of_kp_1: str, rc_of_from: str, rc_of_to: str):
        if kmer_from not in self.nodes:
            self.nodes[kmer_from] = Node(sequence=kmer_from, rc_sequence=rc_of_from)
        if kmer_to not in self.nodes:
            self.nodes[kmer_to] = Node(sequence=kmer_to, rc_sequence=rc_of_to)

        edge = Edge(sequence=kp1_mer_edge, src=self.nodes[kmer_from], dest=self.nodes[kmer_to],
                    kp1=self.k + 1, rev_comp_sequence=rc_of_kp_1)

        self.nodes[kmer_to].path_from[edge] = self.nodes[kmer_from]
        self.nodes[kmer_from].path_to[edge] = self.nodes[kmer_to]

        self.edges[kp1_mer_edge] = edge

    def _fasta_option(self, file):
        with open(file=file) as f_read:
            f_read.readline()
            sequence = []
            for line in tqdm(f_read):
                if line.startswith(DBGraph.FASTA_HEADER):
                    sequence = ''.join(sequence)
                    self._process_kmers_from_sequence(sequence)
                    sequence = []

                else:
                    sequence.append(line.strip())
            else:
                sequence = ''.join(sequence)
                self._process_kmers_from_sequence(sequence)

    def _fastq_option(self, file):

        with open(file=file) as f_read:
            for i, line in tqdm(enumerate(f_read, 1)):
                if i % 4 == 2:
                    self._process_kmers_from_sequence(line.strip())

    def _initialize(self, path_to_fasta, data_format):

        init_option = {"fasta": self._fasta_option, "fastq": self._fastq_option}

        try:
            init_option[data_format](path_to_fasta)

        except KeyError:
            raise FormatError("Provided format isn't supported or incorrect. Check starting command!")

    def _init_get_rev_comp(self, kp1_mer: str):
        kp1_revcomp = self.get_rev_comp(kp1_mer)

        return kp1_revcomp[1:], kp1_revcomp[:-1], kp1_revcomp

    def get_rev_comp(self, sequence):
        reverse_iterator = range(len(sequence) - 1, -1, -1)
        revcomp = ''.join(DBGraph.REVCOMP_DICT[sequence[i]] for i in reverse_iterator)

        self.rev_com_dict[sequence] = revcomp
        self.rev_com_dict[revcomp] = sequence

        return revcomp

    def _simplify(self, node: Node):

        edge_from_prev_node, prev_node = tuple(node.path_from.items())[0]  # as it`s represented like ((Node: Edge), )
        edge_to_next_node, next_node = tuple(node.path_to.items())[0]
        new_edge_sequence = f"{edge_from_prev_node.sequence}{edge_to_next_node.sequence[self.k:]}"

        # appendix_of_edge_to_next_node_coverage = self._compute_edge_coverage(
        #     edge_to_next_node.sequence[self.k:]) // len(edge_to_next_node)

        # new_edge_coverage = (self._compute_edge_coverage(edge_from_prev_node.sequence) // len(
        #     edge_from_prev_node) + appendix_of_edge_to_next_node_coverage) // (
        #                                 len(edge_to_next_node) + edge_from_prev_node)

        self.edges[new_edge_sequence] = Edge(sequence=new_edge_sequence, src=prev_node, dest=next_node, kp1=self.k + 1,
                                             rev_comp_sequence=self.get_rev_comp(new_edge_sequence),
                                             coverage=self._compute_edge_coverage(edge_sequence=new_edge_sequence))

        prev_node.path_to[self.edges[new_edge_sequence]] = next_node
        next_node.path_from[self.edges[new_edge_sequence]] = prev_node

        del prev_node.path_to[edge_from_prev_node]
        del next_node.path_from[edge_to_next_node]
        del self.edges[edge_from_prev_node.sequence]
        del self.edges[edge_to_next_node.sequence]
        del self.nodes[node.sequence]

    def _compute_edge_coverage(self, edge_sequence):
        n = len(edge_sequence)
        coverage = sum(self.frequencies[edge_sequence[i:i + self.k + 1]] + self.frequencies.get(
            self.rev_com_dict[edge_sequence[i:i + self.k + 1]], 0) for i in range(n - self.k))

        return coverage

    def assign_coverage_to_edges(self):
        for edge in tqdm(self.edges):
            self.edges[edge].coverage = self._compute_edge_coverage(edge_sequence=edge)

    def condense_graph(self):

        nodes_to_visit = set(self.nodes.values())
        while nodes_to_visit:
            node = nodes_to_visit.pop()
            if node.is_passing():
                self._simplify(node)

    def remove_alone_node(self, node: Node):
        del self.nodes[node.sequence]

    def remove_edge(self, edge: Edge):
        try:
            edge.src.path_to.pop(edge)
            edge.dest.path_from.pop(edge)

            for node in (edge.dest, edge.src):
                if node.is_alone():
                    self.remove_alone_node(node)
                elif node.is_passing():
                    self._simplify(node)

            del self.edges[edge.sequence]
        except KeyError:
            breakpoint()

    def _remove_inward_tip(self, edge: Edge):

        # breakpoint()
        del edge.dest.path_from[edge]
        del self.nodes[edge.src.sequence]
        del self.edges[edge.sequence]

        if edge.dest.is_passing():
            # breakpoint()
            self._simplify(edge.dest)

        elif edge.dest.is_alone():
            self.remove_alone_node(edge.dest)

    def _remove_outward_tip(self, edge: Edge):

        del edge.src.path_to[edge]
        del edge.dest.path_from[edge]
        del self.nodes[edge.dest.sequence]
        del self.edges[edge.sequence]

        if edge.src.is_passing():
            self._simplify(edge.src)

        elif edge.src.is_alone():
            self.remove_alone_node(edge.src)

    def _is_tip(self, edge: Edge, condition, k_multiplier_for_length=2):
        try:
            if condition(edge) and len(edge.sequence) < k_multiplier_for_length * self.k:

                if all((edge.dest.deg_out == 0, edge.dest.deg_in == 1, edge.src.deg_out > 1)):
                    return 1  # Option 1

                elif all((edge.src.deg_out == 1, edge.src.deg_in == 0, edge.dest.deg_out > 1)):
                    return 2  # Option 2

            else:
                return 0
        except TypeError:
            breakpoint()

    def remove_tips(self, condition, k_multiplier_for_length=2):
        complement_action = {1: 2, 2: 1, 0: 0}
        perform_action = {1: self._remove_outward_tip, 2: self._remove_inward_tip, 0: lambda x: None}

        while True:

            TIPS_FOUND = 0
            visited_edges = set()
            edges = list(self.edges.values())
            for edge in tqdm(edges):
                if edge not in visited_edges and edge.sequence in self.edges:
                    visited_edges.add(edge)

                    action_required = self._is_tip(edge, condition=condition,
                                                   k_multiplier_for_length=k_multiplier_for_length)
                    if action_required:
                        TIPS_FOUND += 1
                        perform_action[action_required](edge)

                    revcomp = self.rev_com_dict[edge.sequence]

                    if action_required and revcomp in self.edges:
                        # complement_action_required = self._is_tip(self.edges[revcomp])
                        # visited_edges.add(self.edges[revcomp])
                        #
                        # if complement_action_required:
                        #     TIPS_FOUND += 1
                        #     perform_action[complement_action_required](self.edges[revcomp])
                        perform_action[complement_action[action_required]](self.edges[revcomp])

            if not TIPS_FOUND:
                break

    def iterative_edges_removing(self, condition):

        visited_edges = set()

        while True:
            BAD_EDGES_FOUND = 0
            edges_items = set(self.edges.items()) - visited_edges

            while edges_items:
                edge_sequence, edge = edges_items.pop()
                revcomp_seq = self.rev_com_dict.get(edge_sequence, None)
                revcomp_edge = self.edges.get(revcomp_seq, None)

                edges_items.discard((revcomp_seq, revcomp_edge))

                visited_edges |= {(edge_sequence, edge), (revcomp_seq, revcomp_edge)}

                if edge_sequence in self.edges and condition(edge):
                    BAD_EDGES_FOUND += 1
                    self.remove_edge(edge)
                    if not edge.is_self_complement and revcomp_edge:
                        self.remove_edge(revcomp_edge)

            if not BAD_EDGES_FOUND:
                break

    def remove_bulges(self):

        self.iterative_edges_removing(condition=self.bad_edge_coverage)

        nodes_directed_connections = defaultdict(list)
        for edge in tqdm(self.edges.values()):
            nodes_directed_connections[(edge.src, edge.dest)].append(edge)

        for bunch_of_connections_between_node_i_node_j in tqdm(nodes_directed_connections.values()):
            for edge in sorted(bunch_of_connections_between_node_i_node_j, key=lambda x: x.coverage)[:-1]:
                self.remove_edge(edge)

        self.remove_tips(condition=lambda x: True, k_multiplier_for_length=4)

        mean_length = np.mean([len(e) for e in self.edges])
        logger.info(f"Mean length {mean_length}")
        self.iterative_edges_removing(condition=lambda edge: all((len(edge.sequence) < 3 * self.k,
                                                                  len(edge.sequence) < mean_length / 100,
                                                                  edge.src not in set(edge.dest.path_to.values()))))

    def write_edges(self):
        saving_file = os.path.join(self.outdir, "edges.fasta")

        logger.info(f"Saving edges to {saving_file}")

        with open(saving_file, "w") as f_write:
            i = 1
            header_template = ">Edge_{}_COV_{}_LEN_{}\n{}\n"
            for edge in tqdm(self.edges.values()):
                f_write.write(header_template.format(i, edge.coverage, len(edge.sequence), edge.sequence))

                i += 1

    def draw_and_save_graph(self):
        dot_file = os.path.join(self.outdir, "graph.dot")
        drawing_file = os.path.join(self.outdir, "graph_planar_view.png")

        logger.info(f"Saving graph structure to {dot_file}")

        G = nx.MultiDiGraph(directed=True)

        # for edge in graph.edges.values():
        #     start, finish, coverage = edge.src, edge.dest, edge.coverage
        #     G.add_edge(start, finish, label=coverage)
        # nx.drawing.nx_pydot.write_dot(G, os.path.join(outdir, "graph.dot"))
        #

        for node in self.nodes.values():
            G.add_node(node, label="")

        for edge in self.edges.values():
            start, finish, coverage = edge.src, edge.dest, edge.coverage
            G.add_edge(start, finish, label=f"COV={coverage} | LEN={len(edge.sequence)}")

        nx.drawing.nx_pydot.write_dot(G, dot_file)

        logger.info(f"Drawing graph to {drawing_file}")
        A = to_agraph(G)
        A.layout("dot")
        A.draw(drawing_file)
