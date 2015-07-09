
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
      index = 0
      Bio::FastaFormat.open(file).each do |entry|
        add_seq index, entry.seq.seq
        index += 1
      end
    end

    def find_starts
      @graph.starts.keys
    end

    # 1) Kmerise contig 1 and build a graph from it. Mark each kmer that
    #    it came from contig 1
    # 2) Kmerise contig 2 and add it to the graph. Find the contig id of all the
    #    contigs that overlap with the kmers from contig 2.
    # 3) Construct a subgraph that contains only the kmers that come from the
    #    contigs found in step 2.
    # 4) Find the start point in this subgraph and traverse it.
    # 5) if there is no fork rename all kmers to match the name of the longest
    #    contig
    # 6) if there is a fork then keep contigs the same
    # 7) um...
    # 8) profit???

    def add_seq id, seq
      kmers = kmerise seq
      # add the kmers as nodes to the graph
      # if they already exist add the contig id to the existing node
      set = Set.new
      kmers.each_with_index do |kmer, index|
        if @graph.exist? (kmer)
          @graph.add_node_value(id, kmer)
          # get the ids of all contigs that these kmers overlap
          contig_list = @graph.get_node_value(kmer)
          contig_list.each do |contig|
            set << contig
          end
        else
          @graph.add([id], kmer)
        end
        # add edge between kmers
        if index > 0
          @graph.add_edge(kmers[index-1], kmer)
        end
      end
      if set.empty?
        # then this contig that was added didn't overlap with anything
        # it was probably the first contig to be added
      else
        # this contig overlapped with another when it was added
        # move to step (3)
        # pull out all the kmers that have the contigs from the set on them
        subgraph = AdjacencyList.new
        @graph.nodes.each do |node_id, list|
          list.value.each do |contig_id|
            if set.include?(contig_id)
              # add graph node to sub graph
              subgraph.add(list.value, node_id)
            end
          end
        end
        # go back through subgraph and get edges from main graph
        # because you can't add edges from nodes that exist to nodes that don't
        # exist yet
        subgraph.nodes.each do |node_id, value|
          neighbour = @graph.edges[node_id]
          if neighbour.size > 0
            neighbour.each do |kmer|
              subgraph.add_edge(node_id, kmer)
            end
          end
        end
        # find start of subgraph
        starts = subgraph.starts
        # hopefully there should be only 1 start
        if starts.size == 1
          # linear = true
          # traverse the subgraph.
          start = starts[0]
          # puts "traversing graph from start point"
          # puts start
          neighbours = subgraph.neighbours(start)
          while neighbours.size > 0 # and linear
            if neighbours.size == 1
              n = neighbours[0]
              neighbours = subgraph.neighbours(n)
            else
              puts "there is a fork in the graph. tines : #{neighbours.size}"
              # linear = false
              # which fork goes the right way
              # the right way has nodes that have `id` on them
              right_way = nil
              neighbours.each do |n|
                list_of_contigs = @graph.get_node_value(n)
                puts "neighbour #{n} has contigs #{list_of_contigs.join(",")}"
                puts "id is #{id}"
                if list_of_contigs.include?(id)
                  #this is the right way
                  right_way = n
                  puts "the right way is #{right_way}"
                  set = set & list_of_contigs
                  puts "set is now #{set.to_a}"
                end
              end
              unless right_way.nil?
                neighbours = subgraph.neighbours(right_way)
              end
              # remove the contigs from the set that were in the other fork
              # of the graph
            end
          end
          # if linear
            # puts "is linear!"
            # rename all the kmers in this subgraph to just come from one
            # contig
            rename = set.to_a.min
            puts "renaming:"
            puts "set is now #{set.to_a}"
            subgraph.nodes.each do |node_id, value|
              # subgraph.set_node_value([rename], node_id)
              # only rename nodes that contain only items in the set
              # and nothing else
              # if set + @graph.get_node_value(node_id) == set
                # puts "setting #{node_id} value from #{@graph.get_node_value(node_id)} to #{rename}"
                # @graph.set_node_value([rename], node_id)
              # end

              # if the node contains all the items in the set
              # remove the items from the set and replace with the min of the
              # set. keep everything else that is already there
              # for example. the node is 1,2,3 and the set is 1,3 then rename
              # is 1 so the node becomes 1,2

              # tmp = list - set.to_a
              # tmp << set.to_a.min

              # if the intersect of the node value and the set contains id
              if (set & @graph.get_node_value(node_id) ).include?(id)
                tmp = @graph.get_node_value(node_id) - set.to_a
                tmp << rename
                tmp.sort!
                puts "setting #{node_id} value from #{@graph.get_node_value(node_id)} to #{tmp}"
                @graph.set_node_value(tmp, node_id)
              end

            end
          # else
            #
          # end
        else
          puts "eek, this shouldn't really happen"
        end

      end
    end

    def extract_seqs
      # puts "extracting all sequences from the graph!"
      # go through all the nodes and get a list of the names of the contigs
      # that are contained within the graph
      set = Set.new

      @graph.nodes.each do |node_id, list|
        list.value.each do |contig_id|
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
        # puts "found first kmer"
        neighbours = @graph.neighbours(start)
        seq = "#{starts[0]}"
        while neighbours.size > 0
          next_kmer = ""
          # puts "there are #{neighbours.size} neighbours"
          neighbours.each do |kmer|
            list = @graph.get_node_value(kmer)
            if list.include?(id)
              next_kmer = kmer
              # puts "#{list.join(",")} includes #{id}. setting next kmer to #{next_kmer}"
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
        # puts ">contig_#{id}"
        # puts seq
        seqs[id]=seq
        # else
          # look through the graph until find a kmer with `id` on it
          # puts "didn't find kmer with contig at the start of the graph"
        # end
      end
      # then traverse the graph for each contig pulling out that sequence
      return seqs
    end

    # __________________________________ ____________________________________

    # def add_seq id, seq
    #   kmers = kmerise seq
    #   (0..kmers.length-1).each do |index|
    #     if index > 0
    #       node_a = kmers[index-1]
    #     end
    #     node_b = kmers[index]
    #     if @graph.exists?(node_b)
    #       puts "adding a value to existing node #{node_b}"
    #       @graph.add_node_value(node_b, id)
    #     else
    #       puts "creating new node #{node_b}"
    #       @graph.add([id], node_b)
    #     end
    #     if index > 0
    #       @graph.add_edge(node_a, node_b)
    #     end
    #   end
    # end


    # def align_seq id, seq
    #   redundant = true
    #   kmers = kmerise seq
    #   first = @graph.get_node_value kmers[0]
    #   hash = {}
    #   if first.size == 0
    #     puts "first size is 0 and so seq is not redundant"
    #     redundant = false
    #     return redundant
    #   else
    #     first.each { |i| hash[i] = 0 }
    #     hash.delete id
    #     # puts "hash with #{id} deleted"
    #     # p hash
    #     kmers.each do |kmer|
    #       puts "counting kmer #{kmer}"
    #       # get the contigs that have this kmer
    #       contigs = @graph.get_node_value(kmer)
    #       # print "contigs: "
    #       # p contigs
    #       # for each contig add 1 to the total counts
    #       contigs.each do |contig|
    #         if hash.key?(contig)
    #           puts "adding 1 to #{contig} in hash"
    #           hash[contig] += 1
    #         end
    #       end
    #       # print "hash: "
    #       # p hash
    #     end
    #     #
    #     hash.each do |k,v|
    #       puts "v: #{v} size: #{kmers.size}"
    #       if v == kmers.size
    #         puts "#{id} is totally contained by #{k}"
    #         return redundant
    #       elsif v < kmers.size
    #         puts "v: #{v} less than kmers: #{kmers.size}"
    #         # have to find out if sequence is like
    #         # ***************------    in this situation you want to merge
    #         # --------*************    0 and 1
    #         #  or:
    #         # *********************    in this situation you want to keep
    #         # ******--------*******    0 and 1 as separate sequences in output
    #       end
    #     end
    #   end
    # end


    def find_path start
      path = []
      puts "start: #{start}"
      path = search path, start
      path
    end

    def search_recursive path, kmer
      puts "path length: #{path.size}\tkmer: #{kmer}"
      neighbours = @graph.adjacent_vertices(kmer)
      if neighbours.size == 0
        path << kmer
        return path
      elsif neighbours.size == 1
        path << kmer
        return search path, neighbours[0]
      elsif neighbours.size > 1
        paths = []
        neighbours.each do |n|
          paths << search(path.dup, n)
        end
        return paths
      end
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
        list << (seq[i..(i+@kmer-1)]).upcase
      end
      # list.each do |s|
      #   puts s
      # end
      list
    end


  end

end