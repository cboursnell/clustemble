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
      @graph.add(0, "ACGT")
      assert_equal 1, @graph.size
    end

    should "add two nodes and connect them with an edge" do
      @graph.add(0, "ACGT")
      @graph.add(0, "CGTT")
      @graph.add_edge("ACGT", "CGTT")
      assert_equal 2, @graph.size, "graph size"
      assert_equal 1, @graph.num_edges, "number of edges"
    end

    should "add two nodes, connect them and get the degree" do
      node1 = "ACGT"
      node2 = "CGTT"
      @graph.add(0, node1)
      @graph.add(0, node2)
      @graph.add_edge(node1, node2)
      assert_equal 1, @graph.in_degree(node2), "in degree"
      assert_equal 1, @graph.out_degree(node1), "out degree"
    end

    should "find start nodes" do
      node0 = "CCGT"
      node1 = "ACGT"
      node2 = "CGTT"
      node3 = "GTTC"
      @graph.add(0, node0)
      @graph.add(0, node1)
      @graph.add(0, node2)
      @graph.add(0, node3)
      @graph.add_edge(node0, node2)
      @graph.add_edge(node1, node2)
      @graph.add_edge(node2, node3)
      starts = @graph.starts
      assert_equal 2, starts.size
      assert starts.include?(node0)
      assert starts.include?(node1)
    end

    should "create nodes with lists as values" do
      @graph.add([], "ACGT")
      assert_equal [], @graph.nodes["ACGT"].value
    end

    should "raise error when adding edges if nodes don't exist" do
      assert_raise RuntimeError do
        @graph.add_edge("ACGT", "CGTG")
      end
    end

    should "test existence of a node" do
      @graph.add(1, "ACGT")
      assert_equal true, @graph.exist?("ACGT")
      assert_equal false, @graph.exist?("cat")
    end

    should "add edge twice and not duplicate" do
      node1 = "ACGT"
      node2 = "CGTT"
      @graph.add(1, node1)
      @graph.add(1, node2)
      @graph.add_edge(node1, node2)
      @graph.add_edge(node1, node2)
      assert_equal 1, @graph.edges[node1].size
    end

    should "get first node with specific id" do
      node1 = "CGTT"
      node2 = "ACGT"
      node3 = "AACG"
      @graph.add([1], node1)
      @graph.add([1], node2)
      @graph.add([1], node3)
      @graph.add_edge(node3, node2)
      @graph.add_edge(node2, node1)
      assert_equal node3, @graph.first_node_with(1), "first node"
    end

  end

end
