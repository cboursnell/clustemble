
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

    def consensus file, name
      add_fasta file
      seqs = extract_seqs
      fasta = ""
      seqs.each.with_index do |(n, seq), i|
        fasta << ">#{name}.#{i+1}\n"
        fasta << "#{seq}\n"
      end
      return fasta
    end

    def add_fasta file
      Bio::FastaFormat.open(file).each do |entry|
        add_seq entry.entry_id, entry.seq.seq
      end
    end

    def find_starts
      @graph.starts.keys
    end

    def add_seq id, seq
      # kmerise the input sequence in an array of kmers
      kmers = kmerise seq
      # puts "kmers size: #{kmers.size}"
      # make a set to store all contigs that `seq` overlaps
      set = Set.new
      kmers.each_with_index do |kmer, index|
        # def add(nodeidentifier, contig, index)
        @graph.add(kmer, id, index)
        # get list of contigs on this kmer
        @graph.get_node_contigs(kmer).each do |contig,index|
          # add contig to set
          set << contig
        end
        # add edge between kmers
        if index > 0
          @graph.add_edge(kmers[index-1], kmer)
        end
      end
      # puts "set: #{set.to_a.join(",")}"
      if set.size == 1 and set & [id] == set
        # this contig that was added doesn't overlap with anything
        # it was probably the first contig to be added
      else
        # get the first kmer that contains a contig from the set on it
        start = @graph.first_node_from_set(set)
        # puts "start: #{start}"
        neighbours = @graph.neighbours(start)
        while neighbours.size > 0
          # print "neighbours "
          # p neighbours
          found = -1
          if neighbours.size == 1
            found = 0
          else
            queue = []
            neighbours.each_with_index do |n, index|
              queue << [n,index]
              # puts "adding #{n} to the queue"
            end
            while queue.size > 0 and found < 0
              front = queue.shift # get the first imte and remove it
              # print "keys: "
              # p @graph.get_node_contigs(front[0]).keys
              if @graph.get_node_contigs(front[0]).keys.include?(id)
                found = front[1]
                # puts "get node contigs includes id #{id}"
              else
                # puts "keep adding to the queue to search deeper"
                nbs = @graph.neighbours(front[0])
                nbs.each do |n|
                  queue << [n, front[1]]
                end
              end
            end
          end
          if found >= 0
            right_way = neighbours[found]
            list_of_contigs = @graph.get_node_contigs(right_way).keys
            if list_of_contigs.include?(id) and neighbours.size > 1
              set = set & list_of_contigs
              # puts "set is now #{set.to_a}"
            end
            neighbours = @graph.neighbours(right_way)
          else
            neighbours = []
            puts "something went wrong"
          end
        end # while neighbours.size > 0
        # puts "set: #{set.to_a}"
        rename = set.to_a.min
        # puts "rename: #{rename}"
        @graph.nodes.each do |node_id, value|
          # puts "kmer #{node_id} contigs: #{@graph.get_node_contigs(node_id).keys}"
          if (set & @graph.get_node_contigs(node_id).keys ).include?(id)
            hash = @graph.get_node_contigs(node_id)
            tmp = hash.keys - set.to_a
            tmp << rename
            # tmp.sort!
            # puts "setting kmer #{node_id} from #{@graph.get_node_contigs(node_id).keys} to #{tmp}"
            new_hash ={}
            tmp.each do |contig_id|
              if hash.key?(contig_id)
                new_hash[contig_id] = hash[contig_id]
              else
                new_hash[contig_id] = 0
              end
            end
            @graph.set_node_contigs(node_id, new_hash)
          end
        end
      end # set.size==1

    end

    def extract_seqs
      # puts "extracting all sequences from the graph!"
      # go through all the nodes and get a list of the names of the contigs
      # that are contained within the graph
      set = Set.new

      @graph.nodes.each do |node_id, node|
        node.contigs.each do |contig_id, index|
          set << contig_id
        end
      end
      # print "set:  "
      # p set
      starts = @graph.starts
      # for each contig id in the graph
      #   find the first kmer that has that id on it
      #   this is not necessarily a kmer with in_degree of 0
      seqs = {}
      set.each do |id|
        # puts "id: #{id}"
        start = @graph.first_node_with id
        # puts "starting at #{start}"
        neighbours = @graph.neighbours(start)
        seq = "#{starts[0]}"
        while neighbours.size > 0
          next_kmer = ""
          # puts "there are #{neighbours.size} neighbours"
          neighbours.each do |kmer|
            hash = @graph.get_node_contigs(kmer)
            if hash.keys.include?(id)
              next_kmer = kmer
              # puts "#{hash.keys.join(",")} includes #{id}. setting next kmer to #{next_kmer}"
            end
          end
          if next_kmer != ""
            seq << next_kmer[-1]
            # puts "seq: #{seq}"
            neighbours = @graph.neighbours(next_kmer)
          else
            neighbours = []
          end
        end
        seqs[id]=seq
      end
      return seqs
    end

    def last_kmer seq
      return seq[-@kmer..-1]
    end

    def kmerise seq
      list = []
      (0..seq.length-@kmer).each do |i|
        kmer = (seq[i..(i+@kmer-1)]).upcase
        list << kmer
      end
      list
    end



  end

end