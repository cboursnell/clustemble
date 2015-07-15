#!/usr/bin/env  ruby

require 'helper'

class TestAdjacencyList < Test::Unit::TestCase

  context 'adjacencylist' do

    setup do
      @graph = Clustemble::AdjacencyList.new
    end

    teardown do
    end

    should "add a node" do
      @graph.add("ACGT", 1, 0)
      assert_equal 1, @graph.size
    end

    should "add two nodes and connect them with an edge" do
      @graph.add("ACGT", 1, 0)
      @graph.add("CGTT", 1, 1)
      @graph.add_edge("ACGT", "CGTT")
      assert_equal 2, @graph.size, "graph size"
      assert_equal 1, @graph.num_edges, "number of edges"
    end

    should "add two nodes, connect them and get the degree" do
      node1 = "ACGT"
      node2 = "CGTT"
      @graph.add(node1, 1, 0)
      @graph.add(node2, 1, 1)
      @graph.add_edge(node1, node2)
      assert_equal 1, @graph.in_degree(node2), "in degree"
      assert_equal 1, @graph.out_degree(node1), "out degree"
    end

    should "add two nodes and then add them again and increment" do
      node1 = "ACGT"
      node2 = "CGTT"
      @graph.add(node1, 1, 0)
      @graph.add(node2, 1, 1)
      @graph.add_edge(node1, node2)
      @graph.add(node1, 2, 0)
      @graph.add(node2, 2, 1)
      assert_equal 2, @graph.nodes[node1].contigs.size
      assert_equal 2, @graph.nodes[node2].contigs.size
    end

    should "find start nodes" do
      node0 = "CCGT"
      node1 = "ACGT"
      node2 = "CGTT"
      node3 = "GTTC"
      @graph.add(node0, 1, 0)
      @graph.add(node1, 1, 1)
      @graph.add(node2, 1, 2)
      @graph.add(node3, 1, 3)
      @graph.add_edge(node0, node2)
      @graph.add_edge(node1, node2)
      @graph.add_edge(node2, node3)
      starts = @graph.starts
      assert_equal 2, starts.size
      assert starts.include?(node0)
      assert starts.include?(node1)
    end

    should "raise error when adding edges if nodes don't exist" do
      assert_raise RuntimeError do
        @graph.add_edge("ACGT", "CGTG")
      end
    end

    should "test existence of a node" do
      @graph.add("ACGT", 1, 0)
      assert_equal true, @graph.exist?("ACGT")
      assert_equal false, @graph.exist?("cat")
    end

    should "add edge twice and not duplicate" do
      node1 = "ACGT"
      node2 = "CGTT"
      @graph.add(node1, 1, 0)
      @graph.add(node2, 1, 1)
      @graph.add_edge(node1, node2)
      @graph.add_edge(node1, node2)
      assert_equal 1, @graph.edges[node1].size
    end

    # should "get first node with specific id" do
    #   node1 = "CGTT"
    #   node2 = "ACGT"
    #   node3 = "AACG"
    #   node4 = "AACG"
    #   @graph.add(Clustemble::Node.new(0, 1, [1]), node1)
    #   @graph.add(Clustemble::Node.new(1, 1, [1]), node2)
    #   @graph.add(Clustemble::Node.new(2, 1, [1]), node3)
    #   @graph.add(Clustemble::Node.new(0, 1, [2]), node4)
    #   @graph.add_edge(node3, node2)
    #   @graph.add_edge(node2, node1)
    #   assert_equal node3, @graph.first_node_with(1), "first node"
    # end

  end

end
