module Clustemble

  # Implements an Adjacency list with indexed nodes
  class AdjacencyList

    ALNode = Struct.new(:value)

    attr_accessor :nodes, :edges

    # Returns a new AdjacencyList
    def initialize
      @nodes = {}
      @edges = {}
      @back_edges = {}
    end

    # Assignment - adds a new node with +:value+, and
    # +:nodeidentifier+, and optionally an array of
    # identifiers of other nodes defining +:edges+.
    # Returns self, so that assignments can be chained.
    def add(value, nodeidentifier)
      node = ALNode.new(value)
      @nodes[nodeidentifier] = node
      @edges[nodeidentifier] = []
      @back_edges[nodeidentifier] = []
      self
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
    def get_node_value nodeidentifier
      @nodes[nodeidentifier].value
    end

    # Set with value of node at +:nodeidentifier+ to +:value+
    def set_node_value(value, nodeidentifier)
      @nodes[nodeidentifier].value = value
    end

    def add_node_value(value, nodeidentifier)
      if @nodes[nodeidentifier].value.is_a?(Array)
        @nodes[nodeidentifier].value << value
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

    def first_node_with id
      first = nil
      @nodes.each do |node_id, list|

        list.value.each do |contig_id|
          if contig_id==id and first.nil?
            first = node_id
          end
        end
      end
      # then check if there are any edges that go from a node to this node
      # and the from node has `id` on it

      # using back edges trace from this node backwards

      found = true
      while found
        found = false
        previous = @back_edges[first]
        if previous.size == 1
          if @nodes[previous[0]].value.include?(id)
            first = previous[0]
            found = true
          end
        else
          previous.each do |n|
            if @nodes[n].value.include?(id)
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
        s += "#{identifier} (#{node.value}) => #{@edges[identifier]} \n"
      end
      s
    end

  end
end