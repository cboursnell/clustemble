
require 'rgl/adjacency'
require 'rgl/bidirectional'

module Clustemble

  class Clustemble

    def initialize kmer
      @nodes = {} # hash
      @graph = RGL::DirectedAdjacencyGraph.new
      @kmer = kmer # size of kmers
    end

    def add_fasta file

    end



  end

end