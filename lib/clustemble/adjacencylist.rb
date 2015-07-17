module Clustemble

  require 'set'

  # Implements an Adjacency list with indexed nodes

  class Node

    attr_accessor :max_count, :visit_count, :contigs, :visit_count

    def initialize contig, index
      # @index = index
      # @max_count = 1   # i think that count might have to be specific to each
      @visit_count = 0 # contig.
      @contigs = Hash.new
      @contigs[contig] = index
    end

    def add_contig(contig, index)
      if @contigs.key?(contig)
      else
        @contigs[contig] = index
        # @max_count += 1
      end
    end

  end

  class AdjacencyList

    attr_accessor :nodes, :edges

    # Returns a new AdjacencyList
    def initialize
      @nodes = {}
      @edges = {}
      @back_edges = {}
    end

    def add(nodeidentifier, contig, index)
      if @nodes.key?(nodeidentifier)
        @nodes[nodeidentifier].add_contig(contig, index)
      else
        @nodes[nodeidentifier] = Node.new(contig, index)
        @edges[nodeidentifier] = []
        @back_edges[nodeidentifier] = []
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

    def get_node_visited nodeidentifier
      if @nodes.key?(nodeidentifier)
        return @nodes[nodeidentifier].visit_count
      else
        return nil
      end
    end

    def get_node_count nodeidentifier
      @nodes[nodeidentifier].max_count
    end

    def inc_node_visited nodeidentifier
      @nodes[nodeidentifier].visit_count += 1
    end

    def reset_node_visited nodeidentifier
      @nodes[nodeidentifier].visit_count = 0
    end

    # Set with value of node at +:nodeidentifier+ to +:value+
    def set_node_contigs(nodeidentifier, contigs)
      if @nodes.key?(nodeidentifier)
        @nodes[nodeidentifier].contigs = contigs
        return true
      else
        return false
      end
    end

    def add_node_contig(id, nodeidentifier)
      if @nodes[nodeidentifier].contigs.is_a?(Array)
        @nodes[nodeidentifier].contigs << id
        # puts "adding #{value} to #{nodeidentifier}"
      else
        raise RuntimeError.new "can't add value to non Array"
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

    def first_node_from_set set
      if !set.is_a?(Set)
        tmpset = Set.new
        tmpset << set
        set = tmpset
      end
      first = nil
      @nodes.each do |node_id, node|
        # puts "#{node_id}\t#{node.contigs.join(",")}"
        node.contigs.each do |contig_id, index|
          # puts "set:#{set.to_a.join(",")}\tcontig: #{contig_id}"
          if set.include?(contig_id) and first.nil?
            first = node_id
          end
        end
      end
      # then check if there are any edges that go from a node to this node
      # and the from node has `id` on it

      # using back edges trace from this node backwards
      abort "shit" if first.nil?
      puts "earliest node found so far is #{first}"
      found = true
      while found
        found = false
        previous = @back_edges[first]
        print "before #{first}. previous: "
        p previous
        node_index = @nodes[first].contigs
        puts "node index of #{first} = "
        p node_index
        previous ||= []
        if previous.size > 0
          previous.each_with_index do |n,i|
            puts "  #{(set & @nodes[n].contigs.keys).to_a}"
            puts "  #{@nodes[n].contigs}"
            subset = set & @nodes[n].contigs.keys
            if subset.size > 0
              subset.each do |s|
                if @nodes[n].contigs[s] < node_index[s]
                  puts "#{@nodes[n].contigs[s]} < #{node_index[s]}"
                  first = n
                  found = true
                end
              end
            end
          end
        end
      end

      return first
    end

    def first_node_with id
      # puts "first node with id = #{id}"
      first = nil
      @nodes.each do |node_id, node|
        # puts "#{node_id}\t#{node.join(",")}"
        node.contigs.each do |contig_id,index|
          if contig_id==id and first.nil?
            first = node_id
          end
        end
      end
      # then check if there are any edges that go from a node to this node
      # and the from node has `id` on it

      # using back edges trace from this node backwards
      abort "shit" if first.nil?
      found = true
      while found
        # puts "looping: #{first}"
        found = false
        previous = @back_edges[first]
        # puts "there are #{previous.length} previous nodes"
        if previous.size == 1
          # print previous
          # p @nodes[previous[0]]
          if @nodes[previous[0]].contigs.keys.include?(id)
            first = previous[0]
            found = true
          end
        else
          previous.each do |n|
            if @nodes[n].contigs.keys.include?(id)
              first = n
              found = true
            end
          end
        end
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