
require 'bio'
require 'set'

module Clustemble

  class Clustemble

    attr_reader :graph

    def initialize kmer
      @nodes = {} # hash
      @graph = AdjacencyList.new
      @kmer = kmer # size of kmers
      @contigs = []
    end

    def add_fasta file
      Bio::FastaFormat.open(file).each do |entry|
        add_seq entry.entry_id, entry.seq.seq
      end
    end

    def add_kmers id, seq
      kmers = kmerise id, seq
      kmers.each_with_index do |kmer, index|
        @graph.add(kmer, id, index)
        # add edge between kmers
        unless index == 0
          @graph.add_edge(kmers[index-1], kmer)
        end
      end
      return kmers
    end

    # 1) a ************
    #    b ----****----
    #
    # 2) a ----****----
    #    b ************
    #
    # 3) a ********----
    #    b ----********
    #
    # 4) a ----********
    #    b ********----
    #
    # 5) a ************
    #    b ****----****
    #
    # 6) a ****----****
    #    b ************
    #
    # b is the sequence we're adding
    #
    def align kmers, a, b
      # get the start kmer for the sequence that is being added
      kmer = @graph.starts[b]
      index_b = 0

      path = []
      kmers.each_with_index do |kmer, index_b|

        # if for each kmer in sequence b i get the contigs for
        # the corresponding kmer in sequence a
        # in indices of the kmers has to increase monotonically
        # so even if there are loops are things there will be a path
        # through the kmer indices that increases monotonically
        # what path finding algorithm is this?
        # you would have to find all paths and then choose the one with the
        # lowest score

        contigs = @graph.get_node_contigs(kmer)
        contigs.each do |i|
          if i[:contig]==a
            path[index_b]||=[]
            path[index_b] << i[:index]
          end
        end
      end
      p path

      return 0
    end

    def remove_contig contig
      @graph.nodes.each do |kmer, node|
        node.remove contig
      end
    end

    def extract_seqs
      return {}
    end

    def kmerise id, seq
      list = []
      (0..seq.length-@kmer).each do |i|
        kmer = (seq[i..(i+@kmer-1)]).upcase
        list << kmer
      end
      list << "end#{id}"
      return list
    end



  end

end