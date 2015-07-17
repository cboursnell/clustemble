
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
      puts "Adding sequence. id : #{id}"
      # set visited counts back to 0 for all nodes
      @graph.nodes.each do |node_id, node|
        node.visit_count = 0
      end
      # kmerise the input sequence in an array of kmers
      kmers = kmerise seq
      puts "kmers size: #{kmers.size}"
      # make a set to store all contigs that `seq` overlaps
      set = Set.new
      subset = Set.new
      kmers.each_with_index do |kmer, index|
        # def add(nodeidentifier, contig, index)
        @graph.add(kmer, id, index)
        # get list of contigs on this kmer
        @graph.get_node_contigs(kmer).each do |contig,index|
          # add contig to set
          set << contig
          subset << contig
        end
        # add edge between kmers
        if index > 0
          puts "edge: #{kmers[index-1]} => #{kmer}"
          @graph.add_edge(kmers[index-1], kmer)
        end
      end
      puts "set: #{set.to_a.sort.join(",")}"
      if set.size == 1 and set & [id] == set
        # this contig that was added doesn't overlap with anything
        # it was probably the first contig to be added
      else
        # get the first kmer that contains a contig from the set on it
        start = @graph.first_node_from_set(set)
        puts "start: #{start}"
        node_index = @graph.get_node_contigs(start)[id]
        puts "start index: #{node_index}"

        neighbours = @graph.neighbours(start)
        while neighbours.size > 0
          print "neighbours: "
          str = []
          neighbours.each do |n|
            str << "#{@graph.get_node_contigs(n)}"
          end
          puts str.join(", ")
          found = -1
          if neighbours.size == 1
            found = 0
          else
            queue = []
            neighbours.each_with_index do |n, index|
              if @graph.get_node_visited(n) == 0
                queue << [n,index]
                print "adding #{n} to the queue "
                print "because #{@graph.get_node_visited(n)} < #{@graph.out_degree(n)}.  "
                print "current: #{@graph.get_node_contigs(start)}\t"
                print "#{@graph.get_node_contigs(n)}"
              elsif @graph.out_degree(n) == 1
                queue << [n,index]
              else
                puts "not adding #{n} because we've already been there"
              end
            end
            while queue.size > 0 and found < 0
              front = queue.shift # get the first item and remove it

              if @graph.get_node_contigs(front[0]).keys.include?(id) and
                 @graph.get_node_contigs(front[0])[id] == node_index+1
                found = front[1]
                puts "found neighbour with id and node_index+1"
              elsif @graph.get_node_contigs(front[0]).keys.include?(id)
                found = front[1]
                puts "found neighbour with id"
              else
                puts "keep adding to the queue to search deeper"
                nbs = @graph.neighbours(front[0])
                nbs.each do |n|
                  queue << [n, front[1]]
                end
              end
            end
          end
          if found >= 0
            start = neighbours[found]
            node_index = @graph.get_node_contigs(start)[id]
            @graph.inc_node_visited(start)
            puts "node #{start} set visited to #{@graph.get_node_visited(start)}"
            list_of_contigs = @graph.get_node_contigs(start).keys
            if list_of_contigs.include?(id) and neighbours.size > 1
              set = set & list_of_contigs
              puts "set is now #{set.to_a}"
            end
            neighbours = @graph.neighbours(start)
          else
            neighbours = []
            puts "something went wrong"
          end
        end # while neighbours.size > 0
        puts "subset: #{subset.to_a}"
        puts "set: #{set.to_a}"
        rename = set.to_a.min
        puts "rename: #{rename}"
        @graph.nodes.each do |node_id, value|
          puts "kmer #{node_id} contigs: #{@graph.get_node_contigs(node_id).keys}"
          if (subset & @graph.get_node_contigs(node_id).keys ).include?(id)
            hash = @graph.get_node_contigs(node_id)
            tmp = hash.keys - subset.to_a
            tmp << rename
            # tmp.sort!
            puts "  setting kmer #{node_id} from #{@graph.get_node_contigs(node_id).keys} to #{tmp}"
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
      puts "extracting all sequences from the graph"
      # go through all the nodes and get a list of the names of the contigs
      # that are contained within the graph
      set = Set.new

      @graph.nodes.each do |node_id, node|
        node.contigs.each do |contig_id, index|
          set << contig_id
        end
      end
      print "set:  "
      p set.to_a
      # starts = @graph.starts
      # for each contig id in the graph
      #   find the first kmer that has that id on it
      #   this is not necessarily a kmer with in_degree of 0
      seqs = {}
      set.each do |id|
        reset # sets node visited counts to 0
        puts "id: #{id}"
        start = @graph.first_node_from_set id
        puts "starting at #{start}"
        neighbours = @graph.neighbours(start)
        node_index = @graph.get_node_contigs(start)[id]
        seq = "#{start}"
        while neighbours.size > 0
          next_kmer = ""
          puts "there are #{neighbours.size} neighbours from #{start}"
          lowest = 1e6
          neighbours.each do |kmer|
            hash = @graph.get_node_contigs(kmer)
            print "  #{kmer} contigs "
            p hash
            puts "visited: #{@graph.get_node_visited(kmer)}. out_degree: #{@graph.out_degree(kmer)}"
            if @graph.out_degree(start) == 1
              lowest = hash[id]
              next_kmer = kmer
            elsif hash.keys.include?(id) and hash[id] < lowest and
              @graph.get_node_visited(kmer) == 0
              lowest = hash[id]
              next_kmer = kmer
              puts "    includes #{id}. setting next kmer to #{next_kmer}"
              puts "    node index found : #{hash[id]}"
            end
          end
          if next_kmer != ""
            seq << next_kmer[-1]
            # puts "seq: #{seq}"
            neighbours = @graph.neighbours(next_kmer)
            @graph.inc_node_visited(next_kmer)
            puts "setting visited count of #{next_kmer} to #{@graph.get_node_visited(next_kmer)}"
            start = next_kmer
            node_index = @graph.get_node_contigs(start)[id]
          else
            neighbours = []
            puts "couldn't find next kmer. setting neighbours to empty list"
          end
        end
        seqs[id] = seq
      end
      return seqs
    end

    def reset
      @graph.nodes.each do |id, node|
        node.visit_count = 0
      end
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