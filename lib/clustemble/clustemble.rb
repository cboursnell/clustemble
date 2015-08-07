
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

    def add_seq id, seq
      @contigs[id] = seq
      if id > 0 # compare with other sequences
        (0..id-1).each do |i|
          puts "comparing #{i} and #{id}"
          kmers1, dups1 = kmerise i, @contigs[i]
          kmers2, dups2 = kmerise id, seq
          if dups1 or dups2
            puts "do more complex alignment"
            sw_align kmers1, kmers2
          else
            puts "do easier alignment"
            sw_align kmers1, kmers2
          end
        end
      end
    end

    def sw_align kmers1, kmers2
      matrix = []
      (0..kmers1.length).each do |x|
        matrix[x] ||= []
        matrix[x][0] = {:score => 0, :arrow => nil}
      end
      (0..kmers2.length).each do |y|
        matrix[0] ||= []
        matrix[0][y] = {:score => 0, :arrow => nil}
      end

      (1..kmers1.length).each do |x|
        (1..kmers2.length).each do |y|
          score_ins = matrix[x][y-1][:score] + 0
          score_match = matrix[x-1][y-1][:score] + 1
          score_del = matrix[x-1][y][:score] + 0
          arrow = nil
          score =-1
          if score_match >= score_ins and score_match >= score_del
            arrow = "M"
            score = score_match
          elsif score_in >= score_match and score_in >= score_del
            arrow = "I"
            score = score_ins
          elsif score_del >= score_match and score_del >= score_ins
            arrow = "D"
            score = score_del
          end
          matrix[x][y] = {:score => score, :arrow => arrow}
        end
      end

      (0..kmers1.length).each do |x|
        (0..kmers2.length).each do |y|
          a = matrix[x][y][:score]
          if a.nil?
            print " - "
          else
            print " #{a} "
          end
        end
        print "\n"
      end

      # find the largest number in the matrix

      max = 0
      startx=0
      starty=0
      (0..kmers1.length).each do |x|
        (0..kmers2.length).each do |y|
          if matrix[x][y][:score] > max
            max = matrix[x][y][:score]
            startx = x
            starty = y
          end
        end
      end
      puts "startx: #{startx}\tstarty: #{starty}\tmax: #{max}"
    end

    def remove_contig contig
      @graph.nodes.each do |kmer, node|
        node.remove contig
      end
    end

    def add_kmers id, seq
      kmers,dups = kmerise id, seq
      set = Set.new
      kmers.each_with_index do |kmer, index|
        @graph.add(kmer, id, index)
        @graph.get_node_contigs(kmer).each do |contig|
          # add contig to set
          set << contig[:contig]
        end
        # add edge between kmers
        if index == 0
          @start_kmer = kmer
        else
          @graph.add_edge(kmers[index-1], kmer)
        end

      end
      return [set,kmers]
    end

    def extract_seqs
      return {}
    end

    def kmerise id, seq
      hash = {}
      list = []
      (0..seq.length-@kmer).each do |i|
        kmer = (seq[i..(i+@kmer-1)]).upcase
        list << kmer
        hash[kmer]||=0
        hash[kmer]+=1
      end
      list << "end#{id}"
      duplicated=false
      hash.each do |k,v|
        if v>1
          duplicated=true
        end
      end
      [list,duplicated]
    end



  end

end