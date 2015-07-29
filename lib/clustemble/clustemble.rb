
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

    def add_seq_old id, seq
      puts "Adding sequence. id : #{id}"
      # set visited counts back to 0 for all nodes
      @graph.nodes.each do |node_id, node|
        node.visit_count = 0
      end
      # kmerise the input sequence in an array of kmers
      kmers = kmerise id, seq
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

    def add_seq id, seq
      verbose = true
      puts "Adding sequence. id : #{id}" if verbose
      kmers = kmerise id, seq
      puts "Kmer count: #{kmers.size}" if verbose
      set = Set.new
      set_all = Set.new
      kmers.each_with_index do |kmer, index|
        @graph.add(kmer, id, index)
        # get list of contigs on this kmer
        @graph.get_node_contigs(kmer).each do |contig|
          # add contig to set
          set << contig[:contig]
          set_all << contig[:contig]
        end
        # add edge between kmers
        if index > 0
          puts "  adding edge #{kmers[index-1]} => #{kmer}" if verbose
          @graph.add_edge(kmers[index-1], kmer)
        end
      end
      puts "set: #{set.to_a.sort.join(",")}" if verbose
      if set.size == 1 and set & [id] == set
        puts "set is size 1. probably adding the first contig" if verbose
      else
        path = [] # list of kmers traversed. used for renaming at the end
        puts "set is size is > 1. overlap this contig with others" if verbose
        start = @graph.first_node set
        path << start
        puts "start is found to be #{start}" if verbose
        info = @graph.get_node_contigs(start)
        # get the lowest index from the info. there should be one with 0
        # but it won't necessarily be from 'id'
        index = 1e6
        current = -1
        others = Hash.new
        info=info.sort_by {|x| [x[:index],-x[:contig]]}
        p info
        info.each do |h|
          print "  each info hash:  " if verbose
          p h if verbose
          i = h[:index]
          puts "  index: #{i}" if verbose
          if i < index
            index = i
            current = h[:contig]
          end
          if i == 0
            others[h[:contig]] = i
          end
        end
        others = others.delete_if {|k,v| k==current}
        puts "current: #{current}" if verbose
        abort "this isn't supposed to happen" if current == -1
        print "others: " if verbose
        p others if verbose
        others_list = []

        # get neighbours
        neighbours = @graph.neighbours(start)
        r = @graph.set_visited(start, current, index, true)
        puts "setting visited on #{start} #{current} #{index}" if verbose
        others.each do |cont, indx|
          puts "   also setting visited on #{cont} #{indx}" if verbose
          r = @graph.set_visited(start, cont, indx, true)
          puts "     output == '#{r}'" if verbose
          if r
            others_list << {:contig => cont, :index => indx, :kmer => start}
          else
            puts "    **** couldn't set visibility on contig: #{cont}, idx: #{indx}"  if verbose
          end
        end
        print "neighbours:" if verbose
        p neighbours if verbose
        count = 0
        while neighbours.size > 0
          prev = start
          candidates = Set.new
          neighbours.each_with_index do |n, n_index|
            n_info = @graph.get_node_contigs(n)
            print "n_info for nbour == #{n_index} : " if verbose
            p n_info if verbose
            n_info.each do |n_hash|
              n_contig_id = n_hash[:contig]
              list_of_contigs << n_contig_id
              n_kmer_id = n_hash[:index]
              if n_contig_id == current and
                 n_kmer_id == index + 1 and
                 !@graph.visited?(n, n_contig_id, n_kmer_id) and
                 n =~ /[ACGTU]+/
                score = 0
              elsif n_contig_id != current and
                    !@graph.visited?(n, n_contig_id, n_kmer_id) and
                    n =~ /[ACGTU]+/

                score = 1
              elsif n =~ /end/
                score = 2
              end
              # if neighbour n doesn't contain id on it somewhere then penalise
              # all candidates of that neighbour

              candidates << { :kmer => n, :index => n_kmer_id,
                              :contig => n_contig_id, :score => score }
              if n_kmer_id == 0 and n_contig_id != current and n_hash[:visit]==false
                unless others.key?(n_contig_id)
                  puts "found a new 'other'" if verbose
                  others[n_contig_id]=n_kmer_id-1
                else
                  puts "somehow found 'new' other that we've seen before" if verbose

                end
              end
            end

          end
          print "list of contigs: "
          p list_of_contigs
          if neighbours.size > 1
            # if there is more than one more neighbour
            # and one of the neighbours has id on it
            # and the other ones don't have id on them
            # then delete the candidates for the neighbours without id
            # OR
            # change the score so that those candidates are penalised
          end

          # if there is a choice between neighbours that have id on them
          # and neighbours that don't, remove the candidates for the neighbours
          # that don't

          print "  candidates: " if verbose
          p candidates if verbose
          if candidates.size >= 2
            r = candidates.to_a.sort_by! { |x| [ x[:score], x[:index] ] }
            print "  r.first :" if verbose
            if r.first[:contig] != current # changing to a different contig
              others = {} # reset others
              info = @graph.get_node_contigs(r.first[:kmer])
              puts "changing to a new contig" if verbose
              print "info:" if verbose
              p info if verbose
              info.each do |d|
                if d[:contig] != r.first[:contig]
                  others[d[:contig]] = d[:index]
                  puts "adding a new 'other'" if verbose
                end
              end
              print "others: " if verbose
              p others if verbose
            end
            p r.first if verbose
            start = r.first[:kmer]
            index = r.first[:index]
            current = r.first[:contig]
          elsif candidates.size > 0
            start = candidates.to_a.first[:kmer]
            index = candidates.to_a.first[:index]
            current = candidates.to_a.first[:contig]
          else
            puts "oh... not sure why this happened"
          end
          # set the visit flag to true on the kmer
          @graph.set_visited(start, current, index, true)
          puts "  setting #{start} contig:#{current} id:#{index} to visited" if verbose
          unless start =~ /end/
            list_of_contigs = []
            @graph.get_node_contigs(start).each do |c|
              list_of_contigs << c[:contig]
            end
            set = set & list_of_contigs
          end
          puts "set is now: #{set.to_a}"
          others.each do |k,v|
            r = @graph.set_visited(start, k, v+1, true)
            puts "    output == '#{r}'" if verbose
            if r
              puts "    also setting visited on #{k} #{v+1}" if verbose
              others[k] = v + 1
              others_list << {:contig => k, :index => v+1, :kmer => start}
            else
              puts "    **** couldn't set visibility on contig: #{k}, idx: #{v+1}" if verbose
              # go back through the previous ones that we'd set on 'other'
              # and unvisit them
              others_list.each do |hash|
                a = @graph.set_visited(hash[:kmer], hash[:contig],
                                       hash[:index], false)
                if a
                  puts "set #{hash[:kmer]} on #{hash[:contig]} back to false" if verbose
                else
                  puts "couldn't unvisit #{hash[:kmer]} on #{hash[:contig]}" if verbose
                end

              end
              others.delete(k)
              puts "deleting contig #{k} from 'others'" if verbose
            end
          end
          set_all << current
          path << start
          # puts "  start: #{start}\tcurrent: #{current}\tindex: #{index}"
          if start == prev
            puts "something went wrong getting a new start kmer" if verbose
            neighbours = []
          else
            neighbours = @graph.neighbours(start)
          end
          info = @graph.get_node_contigs(start)
          count += 1
          if count > 291
            abort "too much clutter"
          end
          print "neighbours:" if verbose
          p neighbours if verbose
        end
        # set_all is a superset of set
        puts "set is now: #{set.to_a.sort.join(", ")}" if verbose
        rename = set.to_a.min + 100
        puts "going to rename all the kmers in the path to #{rename}" if verbose
        # for each kmer in the path rename them and add new kmer indices
        path.each_with_index do |kmer, index|
          contigs = @graph.get_node_contigs(kmer)
          # remove info from contigs where the contig_id is in 'set_all'
          # add new info with contig_id = rename and index = index
          # rename has to be not in set_all
          p contigs if verbose
          contigs.delete_if { |x| set.include?(x[:contig]) }
          p contigs if verbose
          @graph.add_node_contig(kmer, rename, index)
          p contigs if verbose
          puts "__" if verbose
        end
      end
    end #

 def add_seq_2 id, seq
      verbose = true
      puts "Adding sequence. id : #{id}" if verbose
      kmers = kmerise id, seq
      puts "Kmer count: #{kmers.size}" if verbose
      set = Set.new
      set_all = Set.new
      kmers.each_with_index do |kmer, index|
        @graph.add(kmer, id, index)
        # get list of contigs on this kmer
        @graph.get_node_contigs(kmer).each do |contig|
          # add contig to set
          set << contig[:contig]
          set_all << contig[:contig]
        end
        # add edge between kmers
        if index > 0
          puts "  adding edge #{kmers[index-1]} => #{kmer}"
          @graph.add_edge(kmers[index-1], kmer)
        end
      end
      puts "set: #{set.to_a.sort.join(",")}" if verbose
      if set.size == 1 and set & [id] == set
        puts "set is size 1. probably adding the first contig" if verbose
      else
        path = [] # list of kmers traversed. used for renaming at the end
        puts "set is size is > 1. overlap this contig with others" if verbose
        start = @graph.first_node set
        path << start
        puts "start is found to be #{start}" if verbose
        info = @graph.get_node_contigs(start)
        # get the lowest index from the info. there should be one with 0
        # but it won't necessarily be from 'id'
        index = 1e6
        current = -1
        info.each do |h|
          print "  each info hash:  "
          p h
          i = h[:index]
          puts "  index: #{i}"
          if i < index
            index = i
            current = h[:contig]
          end
        end
        puts "current: #{current}"
        abort "this isn't supposed to happen" if current == -1

        # get neighbours
        neighbours = @graph.neighbours(start)
        @graph.set_visited(start, current, index)
        print "neighbours:"
        p neighbours
        count = 0
        while neighbours.size > 0
          prev = start
          candidates = Set.new
          info.each_with_index do |hash, s_index|
            contig_id = hash[:contig]
            kmer_id = hash[:index]
            neighbours.each_with_index do |n,n_index|
              n_info = @graph.get_node_contigs(n)
              print "n_info for neighbour=#{n_index} : "
              p n_info
              n_info.each do |n_hash|
                n_contig_id = n_hash[:contig]
                n_kmer_id = n_hash[:index]
                print "s_index: #{s_index}\t"
                print "id: #{id}\t"
                print "n_contig_id: #{n_contig_id}\t"
                print "contig_id: #{contig_id}\t"
                print "n_kmer_id: #{n_kmer_id}\t"
                print "kmer_id: #{kmer_id}\t"
                print "index: #{index}\t"
                print "visited: #{@graph.visited?(n, n_contig_id, n_kmer_id)}\n"
                if n_contig_id == id and
                   n_contig_id == contig_id and
                   n_kmer_id == kmer_id + 1 and
                   n_kmer_id == index + 1 and
                   n =~ /[ACGT]+/ and  # kmer isn't an end
                   !@graph.visited?(n, n_contig_id, n_kmer_id)
                  candidates << [n, n_kmer_id, n_contig_id, 0]
                elsif n_kmer_id == 0 and # n'bours kmer index is 0
                      n_contig_id == id # n'bours contig id == id of added seq and
                      n =~ /[ACGT]+/ and  # kmer isn't an end
                      !@graph.visited?(n, n_contig_id, n_kmer_id)
                  candidates << [n, n_kmer_id, n_contig_id, 1]
                elsif n_contig_id == current and
                      n_contig_id == contig_id and
                      n_kmer_id == kmer_id + 1 and
                      n_kmer_id == index + 1 and
                      n =~ /[ACGT]+/ and
                      !@graph.visited?(n, n_contig_id, n_kmer_id)
                  candidates << [n, n_kmer_id, n_contig_id, 2]
                elsif contig_id == n_contig_id and
                      n_kmer_id == kmer_id + 1 and
                      n =~ /[ACGT]+/ and  # kmer isn't an end
                      !@graph.visited?(n, n_contig_id, n_kmer_id)
                  candidates << [n, n_kmer_id, n_contig_id, 3]
                elsif contig_id == n_contig_id and
                      n =~ /[ACGT]+/ and  # kmer isn't an end
                      !@graph.visited?(n, n_contig_id, n_kmer_id)
                  candidates << [n, n_kmer_id, n_contig_id, 4]
                elsif n_contig_id == contig_id and
                      n_kmer_id == kmer_id + 1 and
                      n_contig_id == current
                  candidates << [n, n_kmer_id, n_contig_id, 5]
                end
              end
            end
          end # info.each
          print "  candidates: "
          p candidates.to_a
          if candidates.size >= 2
            r = candidates.to_a.sort_by! { |x| x[3] }
            print "  r.first :"
            p r.first
            start = r.first[0]
            index = r.first[1]
            current = r.first[2]
          elsif candidates.size > 0
            start = candidates.to_a.first[0]
            index = candidates.to_a.first[1]
            current = candidates.to_a.first[2]
          else
            puts "oh... not sure why this happened"
          end
          # set the visit flag to true on the kmer
          @graph.set_visited(start, current, index)
          puts "  setting #{start} #{current} #{index} to visited"
          set_all << current
          path << start
          # puts "  start: #{start}\tcurrent: #{current}\tindex: #{index}"
          if start == prev
            puts "something went wrong getting a new start kmer"
            neighbours = []
          else
            neighbours = @graph.neighbours(start)
          end
          info = @graph.get_node_contigs(start)
          count += 1
          if count > 191
            abort "too much clutter"
          end
          print "neighbours:"
          p neighbours
        end
        # set_all is a superset of set
        puts "set_all is now: #{set_all.to_a.sort.join(", ")}"
        rename = set_all.to_a.min + 100
        puts "going to rename all the kmers in the path to #{rename}"
        # for each kmer in the path rename them and add new kmer indices
        path.each_with_index do |kmer, index|
          contigs = @graph.get_node_contigs(kmer)
          # remove info from contigs where the contig_id is in 'set_all'
          # add new info with contig_id = rename and index = index
          # rename has to be not in set_all
          p contigs
          contigs.delete_if { |x| set_all.include?(x[:contig]) }
          p contigs
          @graph.add_node_contig(kmer, rename, index)
          p contigs
          puts "__"
        end
      end
    end #


    def extract_seqs_old
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
      abort "done"
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

    def extract_seqs
      verbose = true
      seqs = {}
      puts "extracting all sequences from the graph" if verbose
      # go through all the nodes and get a list of the names of the contigs
      # that are contained within the graph
      set = Set.new

      @graph.nodes.each do |kmer, node|
        unless kmer =~ /end/
          node.contigs.each do |hash|
            set << hash[:contig]
            if hash[:contig] < 100
              puts "id of #{hash[:contig]} came from #{kmer}" if verbose
            end
          end
        end
      end
      print "set:  " if verbose
      p set.to_a if verbose
      # for each id in the set
      set.each do |seq_id|
        puts "tracing sequence for #{seq_id}" if verbose
        # go through the graph and find the first node
        start = ""
        current_index = 0
        seq = ""
        @graph.nodes.each do |kmer, node|
          node.contigs.each do |hash|
            contig_id = hash[:contig]
            index = hash[:index]
            if contig_id == seq_id and index == current_index
              start = kmer
            end
          end
        end
        seq << start
        puts "start is #{start}" if verbose
        # start from start and just follow the increasing indices until there
        # are no possible neighbours
        neighbours = @graph.neighbours(start)
        while neighbours.size > 0
          found = false
          # puts "there are #{neighbours.size} neighbours"
          neighbours.each do |kmer|
            unless kmer =~ /end/
              puts "  #{kmer}\t#{@graph.get_node_contigs(kmer)}" if verbose
              @graph.get_node_contigs(kmer).each do |hash|
                if hash[:index] == current_index + 1
                  start = kmer
                  found = true
                  puts "    found next kmer is #{kmer} with index: #{current_index+1}" if verbose
                  current_index += 1
                end
              end
            end
          end
          if found
            seq << start[-1]
            neighbours = @graph.neighbours(start)
          else
            neighbours = []
          end
        end
        seqs[seq_id]=seq
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

    def kmerise id, seq
      list = []
      (0..seq.length-@kmer).each do |i|
        kmer = (seq[i..(i+@kmer-1)]).upcase
        list << kmer
      end
      list << "end#{id}"
      list
    end



  end

end