
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
      verbose = true
      puts "Adding sequence. id : #{id}" if verbose

      # add sequence.
      # - kmerise
      # - create graph
      # - add edges to graph
      set, kmers = add_kmers(id, seq)
      # if the set is just 'id'
      if set.size == 1 and set & [id] == set
        puts "set is size 1. adding the first contig" if verbose
      else
        # start from index 0 of contig id
        # puts "there are #{kmers.size} kmers in seq #{id}"

        # for each kmer check what other contigs are on that kmer
        # it should only be contigs with lower ids
        # there are 3 possible cases
        # 1a) totally redundant with a smaller than b
        # 1b) totally redundant with a larger than b
        # 2a) overlapping
        # 2b) overlapping in two places
        # for example:
        # if you're adding seq1 and there are 10 kmers and all 10 kmers have
        #   seq 0 on them then it is case 1 and the result should be to rename
        #   all the kmers from seq1 to be seq0
        # if you're adding seq1 and there are 10 kmers and 6 kmers have seq0
        #   on them then build up a picture of how they lie
        #   2a) would look like this [*]
        #   2b) would look like this [*-*]

        contigs = {}
        set.to_a.each do |c|
          contigs[c] = []
        end
        contigs.delete(id)
        other_index=-1
        this_index=-1
        # for each of the kmers we've added to the graph...
        kmers.each_with_index do |kmer, index|
          puts "kmer: #{kmer} #{index}"
          unless kmer=~/end/
            # get the info about that kmer in the graph
            list = @graph.get_node_contigs(kmer)
            # for all the other kmers that this sequence overlaps
            print "  "
            p list
            contigs.each do |c, rep|
              unless c==id
                found = false
                other_index = -1
                list.each do |info|
                  # print "  info: "
                  # p info
                  if info[:contig]==c
                    found = true
                    other_index = info[:index]
                  end
                end
                if other_index>=0
                  puts "  difference: #{other_index - index}"
                end
                if found
                  rep << "*"
                else
                  rep << "-"
                end
              end
            end
          end
        end
        p contigs
        cases = {}
        contigs.each do |c, rep|
          exons = 0
          gaps = 0
          s = ""
          rep.each do |i|
            if i!=s
              s = i
              if i!="-"
                exons+=1
              elsif i=="-"
                gaps+=1
              end
            end
          end
          puts "gaps: #{gaps}\texons:#{exons}"
          if (gaps==0 and exons>0) or (gaps==2 and exons==1)
            # puts "  case1"
            cases[c]=1
            # find if c or id are longer

          elsif gaps==1 and exons==1
            # puts "  case2"
            cases[c]=2
          elsif exons>gaps
            # puts "  case3"
            cases[c]=3
          else
            puts "what's case4???"
          end
        end
        # print "cases: "
        # p cases
        mark_this_contig_for_removal = false
        # first go through all cases and find if there are any of type 1
        #   if there are then remove one of the contigs
        # if there are no type 1 then see if there are any type 2
        #   deal with these
        # if theere are types 3 just ignore this
        no_type_one = true
        cases.each do |c, type|
          if type==1
            puts "contig #{c} is type #{type}"
            no_type_one = false
            # go through all kmers that are from contig c and rename them to
            # in order to merge c and id
            # if contig c is smaller remove it
            # if contig id is smaller then remove it
            puts "  contig #{c} is #{self.length(c)} long"
            puts "  contig #{id} is #{self.length(id)} long"
            if self.length(c) < self.length(id)
              # remove c
              puts "  removing contig #{c} because it's shorter than #{id}"
              remove_contig(c)
            else
              # remove id - maybe leave this until after in case there's
              # case 2 later.
              mark_this_contig_for_removal = true
            end
          end
        end

        cases.each do |c, type|
          if type==2 and no_type_one
            puts "contig #{c} is type #{type}"
            # go through all kmers that are from contig c and rename them
            # to
            puts "  renaming contig #{c} to contig #{id}"
            rename(c, id)
          elsif type==3
            # don't do any renaming
          end
        end
        if mark_this_contig_for_removal
          puts "removing this contig #{id}"
          remove_contig id
        end
        # if any of the contigs that id aligned with had a case1
        # then just remove the contig as it is completely redundant


        # from the end trace onwards
        # then from the first kmer trace backwards
        # if the contig is redundant rename all kmers to
      end

    end

    def remove_contig contig
      @graph.nodes.each do |kmer, node|
        node.remove contig
      end
    end

    def rename contig1, contig2
      # used for case2 situations
      # find out which contig is first
      # find a kmer that has both contig 1 and contig 2 on it
      # if the index of contig1 is higher then contig1 is first
      # otherwise contig2 is first
      #
      # the difference should be the same for every overlapping kmer
      # change the contig id on contig1 to contig2
      # then if contig1 is first add then difference to contig2 kmers
      #
      # then if contig2 is first add then difference to contig1 kmers
      #

      # 1) find a kmer that has both contig1 and contig2 on it
      difference = 0
      found_kmer = ""
      @graph.nodes.each do |kmer, node|
        found1=false
        found2=false
        index1=0
        index2=0
        node.contigs.each do |item|
          if item[:contig]==contig1
            found1=true
            index1=item[:index]
          end
          if item[:contig]==contig2
            found2=true
            index2=item[:index]
          end
        end
        if found1 and found2
          puts "#{kmer} has #{contig1} and #{contig2}"
          found_kmer = kmer
          difference = index1-index2
        end
      end
      puts "difference is #{difference}"
      if difference < 0
        puts "contig2 is first"
        # 2) on all kmers that had contig1 and contig2 on them remove contig1
        #      so there isn't any duplication
        @graph.nodes.each do |kmer, node|
          if node.has_all_contigs?([contig1,contig2])
            node.remove contig1
          end
        end
        # on all the kmers that have contig1 on them, increase the index by
        #   difference
        @graph.nodes.each do |kmer, node|
          if node.has_contig?(contig1)
            node.add_to_index(contig1, difference.abs)
          end
        end
        # on all the kmers that are just contig1 rename the contig id to be
        #   contig2
        @graph.nodes.each do |kmer, node|
          if node.has_contig?(contig1) and !node.has_contig?(contig2)
            node.rename(contig1, contig2)
          end
        end
      elsif difference > 0
        puts "contig1 is first"
        # 2) on all kmers that had contig1 and contig2 on them remove contig2
        #      so there isn't any duplication
        @graph.nodes.each do |kmer, node|
          if node.has_all_contigs?([contig1,contig2])
            node.remove contig2
          end
        end
        # on all the kmers that have contig2 on them, increase the index by
        #   difference
        @graph.nodes.each do |kmer, node|
          if node.has_contig?(contig2)
            node.add_to_index(contig2, difference.abs)
          end
        end
        # on all the kmers that are just contig1 rename the contig id to be
        #   contig2
        @graph.nodes.each do |kmer, node|
          if node.has_contig?(contig1) and !node.has_contig?(contig2)
            node.rename(contig1, contig2)
          end
        end

      else
        puts "the contigs are the same"
      end

    end

    def add_kmers id, seq
      kmers = kmerise id, seq
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
            # if hash[:contig] < 100
              # puts "id of #{hash[:contig]} came from #{kmer}" if verbose
            # end
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
                if hash[:index] == current_index + 1 and hash[:contig]==seq_id
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

    def length contig
      count =0
      @graph.nodes.each do |id, node|
        if node.has_contig?(contig)
          count+=1
        end
      end
      return count
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