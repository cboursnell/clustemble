module Clustemble

  require 'set'

  # Implements an Adjacency list with indexed nodes

  class Node

    attr_accessor :kmer, :contigs

    def initialize kmer
      @kmer = kmer
      @contigs = []
    end

    def add_contig(contig, index)
      @contigs <<  { :contig => contig, :index => index, :visit => false }
    end

    def has_contig? contig
      found = false
      @contigs.each do |h|
        if h[:contig]==contig
          found = true
        end
      end
      return found
    end

    def has_contigs? set
      found = false
      @contigs.each do |h|
        if set.include?(h[:contig])
          found = true
        end
      end
      return found
    end

    def indices contig
      list = []
      @contigs.each do |h|
        if h[:contig]==contig
          list << h[:index]
        end
      end
      return list
    end

    def visited? contig, index
      v = false
      @contigs.each do |c|
        if c[:contig]==contig and c[:index]==index
          v = c[:visit]
        end
      end
      return v
    end

    def set_visited contig, index, bool
      set = false
      @contigs.each do |c|
        if c[:contig]==contig and c[:index]==index
          # puts "+ visible set on #{contig}\t#{index}\t#{@kmer}"
          c[:visit] = bool
          set = true
        end
      end
      return set
    end

  end # class Node

  class AdjacencyList

    attr_accessor :nodes, :edges

    # Returns a new AdjacencyList
    def initialize
      @nodes = {}
      @edges = {}
      @back_edges = {}
    end

    def add(kmer, contig, index)
      if @nodes.key?(kmer)
        @nodes[kmer].add_contig(contig, index)
      else
        node = Node.new(kmer)
        node.add_contig(contig, index)
        @nodes[kmer] = node
        @edges[kmer] = []
        @back_edges[kmer] = []
      end
    end

    def add_node(node, nodeidentifier)
      if node.is_a?(Node)
        @nodes[nodeidentifier] = node
        @edges[nodeidentifier] = []
        @back_edges[nodeidentifier] = []
      else
        raise RuntimeError.new "must add Node object"
      end
    end

    # Removal - deletes the node at +:nodeidentifier+, which should be
    # an integer index if this is an indexed adjacency list, or the name
    # of the node if this is a names adjacency list.
    def delete(nodeidentifier)
      node = @nodes[nodeidentifier]
      @nodes[nodeidentifier] = nil
      @edges.delete node
      @back_edges.delete node
    end

    # Removal - deletes the edge(s) +:edges+ connected to the node
    # referenced by +:nodeidentifier+.
    def delete_edge(nodeidentifier, *edges)
      alledges = @edges[nodeidentifier]
      edges.each { |edge| alledges.delete edge }
    end

    # Returns the value of the node with +:nodeidentifier+
    def get_node_contigs nodeidentifier
      if @nodes.key?(nodeidentifier)
        return @nodes[nodeidentifier].contigs
      else
        return nil
      end
    end

    def get_node_index nodeidentifier
      if @nodes.key?(nodeidentifier)
        return @nodes[nodeidentifier].index
      else
        return nil
      end
    end

    def visited? kmer, contig, index
      if @nodes.key?(kmer)
        return @nodes[kmer].visited?(contig, index)
      else
        return nil
      end
    end

    def set_visited kmer, contig, index, bool
      found=false
      if @nodes.key?(kmer)
        found = @nodes[kmer].set_visited contig, index, bool
      end
      return found
    end

    # Set with value of node at +:nodeidentifier+ to +:value+
    def add_node_contig(kmer, id, index)
      if @nodes.key?(kmer)
        @nodes[kmer].add_contig(id, index)
        return true
      else
        return false
      end
    end

    # Adds an edge from node with identifier +:x+ to node
    # with identifier +:y+.
    def add_edge(x, y) # from x to y
      if @nodes.key?(x) and @nodes.key?(y)
        if @edges.key?(x)
          unless @edges[x].include?(y)
            @edges[x] << y
          end
        end
        if @back_edges.key?(y)
          unless @back_edges[y].include?(x)
            @back_edges[y] << x
          end
        end
      else
        raise RuntimeError.new "#{x} and #{y} not both present"
      end
    end

    # True if +:x+ and +:y+ are connected by an edge.
    def adjacent?(x, y)
      @edges[x].include?(y) || @edges[y].include?(x)
    end

    # Return an array of identifiers of all nodes connected to
    # node at +:nodeidentifier+ by edges.
    def neighbours nodeidentifier
      @edges[nodeidentifier]
    end

    def out_degree nodeidentifier
      @edges[nodeidentifier].size
    end

    def in_degree nodeidentifier
      degree=0
      @edges.each do |fromnode, list|
        list.each do |tonode|
          if tonode==nodeidentifier
            degree+=1
          end
        end
      end
      degree
    end

    def starts # return nodes with an in degree of 0
      degrees={}
      @nodes.each do |k,v|
        degrees[k]=0
      end
      @edges.each do |fromnode, list|
        list.each do |tonode|
          degrees[tonode]+=1
        end
      end
      degrees.delete_if {|k,v| v > 0}
      degrees.keys
    end

    def first_node_old set
      # if set is an integer then make it a set
      if !set.is_a?(Set)
        tmpset = Set.new
        tmpset << set
        set = tmpset
      end
      first = ""
      # scan through the nodes and just get the first node it finds
      # p set.to_a
      @nodes.each do |kmer, node|
        if node.has_contigs?(set)
          first = kmer
          break
        end
      end
      # puts "first: #{first}"
      if first==""
        abort "error type 234"
      end
      # now from this kmer work backwards to find the first kmer
      if @back_edges.key?(first)
        previous = @back_edges[first]
        # print "previous: "
        # p previous
        count=0
        while previous.size > 0
          contigs = @nodes[first].contigs
          print "contig info: "
          p contigs
          if previous.size == 1
            print "info on only previous node:  "
            p previous[0]
            p @nodes[previous[0]].contigs
            if @nodes[previous[0]].has_contigs? set
              first = previous[0]
              previous = @back_edges[first]
              # puts "setting first to #{first}"
              # puts "previous size #{previous.size}"
            else
              @nodes[previous[0]].contigs.each do |hash|
                set << hash[:contig]
                puts "adding #{hash[:contig]} to set"
              end
              previous = @back_edges[first]
              puts "set is now: #{set.to_a}"
            end
          else
            puts "looking through #{previous.size} previous nodes"
            list = []
            previous.each do |prev|
              if @nodes[prev].has_contigs?(set)
                list << prev
              end
            end
            puts "made a list of size #{list.size}"
            current = 1e6
            contigs.each do |item|
              print "item1:"
              p item
              i = item[:index]
              if i < current
                current = i
              end
            end
            next_kmer = ""
            puts "set index to #{current}"
            list.each do |kmer|
              @nodes[kmer].contigs.each do |item|
                i = item[:index]
                print "item:"
                p item
                if i == current - 1
                  puts "found a kmer that's before"
                  next_kmer = kmer
                end
              end
            end
            if next_kmer == ""
              next_kmer = list.shuffle.first
              puts "picking a random one... this will work right?"
            end
            first = next_kmer
            previous = @back_edges[first]
          end
          count+=1
          abort "too many times" if count>300
        end # while
      end
      return first
    end

    def first_node set
      # if set is an integer then make it a set
      if !set.is_a?(Set)
        tmpset = Set.new
        tmpset << set
        set = tmpset
      end
      first = ""
      # scan through the nodes and just get the first node it finds
      # p set.to_a
      @nodes.each do |kmer, node|
        if node.has_contigs?(set)
          first = kmer
          break
        end
      end
      # puts "first: #{first}"
      if first==""
        abort "error type 234"
      end
      # count=0
      # now from this kmer work backwards to find the first kmer
      while @back_edges[first].size > 0
        info = @nodes[first].contigs
        index = 1e6
        current = -1
        info.each do |hash|
          if hash[:index] < index
            index = hash[:index]
            current = hash[:contig]
          end
        end
        # puts "first: #{first}"
        previous = @back_edges[first]
        if previous.size == 1
          first = previous[0]
          # puts "setting first to #{first} 1"
        elsif previous.size >= 2
          candidates = []
          # print "previous: "
          # p previous
          previous.each do |prev|
            contigs = @nodes[prev].contigs
            contigs.each do |item|
              c = item[:contig]
              i = item[:index]
              if c == current and index - 1 == i
                candidates << {:kmer=>prev,:contig=>c,:index=>i,:score=>0}
              else
                candidates << {:kmer=>prev,:contig=>c,:index=>i,:score=>1}
              end
            end
          end
          # print "candidates: "
          # p candidates
          res = candidates.sort_by! { |x| [ x[:score], x[:index] ] }
          first = res.first[:kmer]
          # puts "setting first to #{first} 2"
        else
          puts "first: #{first}"
          abort "this shouldn't happen"
        end
        # count+=1
        # abort "schtop!" if count > 100
      end
      return first
    end

    def size
      @nodes.size
    end

    def num_edges
      count=0
      @edges.each do |k,v|
        if v.length > 0
          count+=1
        end
      end
      count
    end

    def exist? nodeidentifier
      @nodes.key?(nodeidentifier)
    end

    # Return a string representation of the graph
    def to_s
      s = ""
      @nodes.each do |identifier, node|
        s += "#{identifier} (#{node.contigs.join(",")}) => #{@edges[identifier]} \n"
      end
      s
    end

  end
end