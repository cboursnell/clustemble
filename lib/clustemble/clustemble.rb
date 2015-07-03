
require 'bio'
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
      Bio::FastaFormat.open(file).each do |entry|
        list = kmerise entry.seq.seq
        list.each_with_index do |kmer, index|
          if index > 0
            @graph.add_edge(list[index-1], list[index])
          end
        end
      end
    end

    def find_start
      start = ""
      in_degree = {}
      @graph.each_vertex do |vertex|
        in_degree[vertex] ||= 0
        @graph.adjacent_vertices(vertex).each do |neighbour|
          in_degree[neighbour] ||= 0
          in_degree[neighbour] += 1
        end
      end
      in_degree.each do |vertex,count|
        if count == 0
          start = vertex
        end
      end
      return start
    end

    def traverse start
      seq = []
      seq << start
      last = false
      while !finished(seq)
        (0..seq.size-1).each do |i|
          kmer = last_kmer seq[i]
          if @graph.has_vertex?(kmer)
            neighbours = @graph.adjacent_vertices(kmer)
            if neighbours.size > 1
              puts "seq #{i}\tkmer #{kmer}\tneighbours: #{neighbours.size}"
              str = seq[i].dup
              # add the last character of the first neighbour
              seq[i] << neighbours[0][-1]
              # add the last character of the next neighbours
              (1..neighbours.size-1).each do |j|
                # and make a new sequence
                seq << "#{str}#{neighbours[j][-1]}"
              end
            elsif neighbours.size == 1
              seq[i] << neighbours[0][-1]
            else
              # puts "last kmer: #{kmer}"
              seq[i] << "|"
            end
          end
        end
      end
      seq
    end

    def finished seq
      done = true
      seq.each do |s|
        if s[-1] != "|"
          done = false
        end
      end
      done
    end

    def last_kmer seq
      return seq[-@kmer..-1]
    end

    # returns a list of kmers
    def kmerise seq
      list = []
      (0..seq.length-@kmer).each do |i|
        list << seq[i..(i+@kmer-1)]
      end
      # list.each do |s|
      #   puts s
      # end
      list
    end


  end

end