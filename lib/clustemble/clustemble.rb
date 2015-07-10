
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
      puts "\nADDING SEQUENCE: ID #{id}"
      kmers = kmerise seq
      puts "kmerised sequence to #{kmers.size} kmers"
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
        puts "first contig added"
        # then this contig that was added didn't overlap with anything
        # it was probably the first contig to be added
      else
        # this contig overlapped with another when it was added
        # move to step (3)
        # pull out all the kmers that have the contigs from the set on them
        subgraph = AdjacencyList.new
        puts "building subgraph. looking for nodes with id: #{id}"
        puts "set is now #{set.to_a}"
        @graph.nodes.each do |node_id, list|
          list.value.each do |contig_id|
            if set.include?(contig_id)
              # add graph node to sub graph
              subgraph.add(list.value, node_id)
            end
          end
        end
        puts "added nodes to subgraph. subgraph contains #{subgraph.size} nodes"
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
          puts "traversing graph from start point: #{start} contigs:#{@graph.get_node_value(start).join(",")}"
          # puts start
          neighbours = subgraph.neighbours(start)
          while neighbours.size > 0 # and linear
            if neighbours.size == 1
              n = neighbours[0]
              neighbours = subgraph.neighbours(n)
              puts "next node is #{n} contigs:#{@graph.get_node_value(n).join(",")}"
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
              else
                puts "couldn't find id:#{id} on immediate neighbours"
                # search forward along both paths of the fork until either
                # 1) the two arms of the fork meet up. (don't try this)
                # 2) a kmer is found with `id` on it (try this!)

                # can probably replace the examination of the neighbours above
                # with this approach below as it basically does the same thing
                queue = []
                # add neighbours to queue with markers to say which is which
                neighbours.each_with_index do |n, index|
                  queue << [n, index]
                end
                found = -1
                while queue.size > 0 and found < 0
                  front = queue.shift # get the first item and remove it
                  puts "   searching: #{front[0]} from path #{front[1]}. contigs #{subgraph.get_node_value(front[0]).join(",")}"
                  if subgraph.get_node_value(front[0]).include?(id)
                    # found it
                    found = front[1]
                    puts "found a kmer with #{id} on it. found=#{found}"
                  else
                    nbs = subgraph.neighbours(front[0])
                    nbs.each do |n|
                      queue << [n, front[1]]
                    end
                  end
                end
                if found >= 0
                  right_way = neighbours[found]
                  neighbours = subgraph.neighbours(right_way)
                else
                  puts "something went really wrong"
                  abort "sorry"
                end
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
            # if the node contains all the items in the set
            # remove the items from the set and replace with the min of the
            # set. keep everything else that is already there
            # for example. the node is 1,2,3 and the set is 1,3 then rename
            # is 1 so the node becomes 1,2

            if (set & @graph.get_node_value(node_id) ).include?(id)
              tmp = @graph.get_node_value(node_id) - set.to_a
              tmp << rename
              tmp.sort!
              # puts "setting #{node_id} value from #{@graph.get_node_value(node_id)} to #{tmp}"
              @graph.set_node_value(tmp, node_id)
            end

          end
          # else
            #
          # end
        else
          puts "eek, this shouldn't really happen"
        end

      end # set.empty?
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



    # def find_path start
    #   path = []
    #   puts "start: #{start}"
    #   path = search path, start
    #   path
    # end

    # def search_recursive path, kmer
    #   puts "path length: #{path.size}\tkmer: #{kmer}"
    #   neighbours = @graph.adjacent_vertices(kmer)
    #   if neighbours.size == 0
    #     path << kmer
    #     return path
    #   elsif neighbours.size == 1
    #     path << kmer
    #     return search path, neighbours[0]
    #   elsif neighbours.size > 1
    #     paths = []
    #     neighbours.each do |n|
    #       paths << search(path.dup, n)
    #     end
    #     return paths
    #   end
    # end

    # def traverse start
    #   seq = []
    #   seq << start
    #   last = false
    #   while !finished(seq)
    #     (0..seq.size-1).each do |i|
    #       kmer = last_kmer seq[i]
    #       if @graph.has_vertex?(kmer)
    #         neighbours = @graph.adjacent_vertices(kmer)
    #         if neighbours.size > 1
    #           puts "seq #{i}\tkmer #{kmer}\tneighbours: #{neighbours.size}"
    #           str = seq[i].dup
    #           # add the last character of the first neighbour
    #           seq[i] << neighbours[0][-1]
    #           # add the last character of the next neighbours
    #           (1..neighbours.size-1).each do |j|
    #             # and make a new sequence
    #             seq << "#{str}#{neighbours[j][-1]}"
    #           end
    #         elsif neighbours.size == 1
    #           seq[i] << neighbours[0][-1]
    #         else
    #           # puts "last kmer: #{kmer}"
    #           seq[i] << "|"
    #         end
    #       end
    #     end
    #   end
    #   seq
    # end

    # def finished seq
    #   done = true
    #   seq.each do |s|
    #     if s[-1] != "|"
    #       done = false
    #     end
    #   end
    #   done
    # end

    def last_kmer seq
      return seq[-@kmer..-1]
    end

    # returns a list of kmers
    def kmerise seq
      hash = {}
      list = []
      (0..seq.length-@kmer).each do |i|
        kmer = (seq[i..(i+@kmer-1)]).upcase
        list << kmer
        if hash.key?(kmer)
          puts "this is bad. the same kmer is in this sequence twice"
          puts "i think this might lead to problems"
        end
        hash[kmer]||=0
        hash[kmer]+=1
      end
      # list.each do |s|
      #   puts s
      # end
      list
    end


  end

end