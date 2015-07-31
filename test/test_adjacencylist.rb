#!/usr/bin/env  ruby

require 'helper'

class TestAdjacencyList < Test::Unit::TestCase

  context 'adjacencylist' do

    setup do
      @graph = Clustemble::AdjacencyList.new
    end

    teardown do
    end

    should "create a node" do
      node = Clustemble::Node.new("ACGT")
      node.add_contig(1, 0)
      assert node.has_contig?(1), "has contig"
      assert_equal [0], node.indices(1), "node indices"
      node.add_contig(1, 33)
      assert_equal [0,33], node.indices(1), "node indices2"
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

    should "get first node with specific id" do
      # seq = AACGTT

      seq1 = "AACGTT"
      k_size = 4
      [seq1].each_with_index do |seq, id|
        kmers = []
        (0..seq.length-k_size).each do |i|
          kmer = (seq[i..(i+k_size-1)]).upcase
          kmers << kmer
        end
        kmers.each_with_index do |kmer, index|
          @graph.add(kmer, id, index)
          @graph.add_edge(kmers[index-1], kmer) if index > 0
        end
      end

      assert_equal "AACG", @graph.first_node(0), "first node"
    end

    should "get the first node in a more complex sequence" do
      seq1 = "TTGGAATCGGTGACCGGCATGAATTTGACAGATT"
      seq2 = "AATTGGAATCGGTGACCGGCATGAATTTGACAGATTGGAATCGGTGACCGGCATGAATTTGACAGTA"
      k_size = 31
      [seq1, seq2].each_with_index do |seq, id|
        kmers = []
        (0..seq.length-k_size).each do |i|
          kmer = (seq[i..(i+k_size-1)]).upcase
          kmers << kmer
        end
        kmers.each_with_index do |kmer, index|
          @graph.add(kmer, id, index)
          @graph.add_edge(kmers[index-1], kmer) if index > 0
        end
      end

      start = @graph.first_node 1
      assert_equal "AATTGGAATCGGTGACCGGCATGAATTTGAC", start

    end

    should "get the first node in a very simple sequence!" do
      seq1 ="TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      seq2 = seq1[4..41]

      k_size = 31
      [seq1, seq2].each_with_index do |seq, id|
        kmers = []
        (0..seq.length-k_size).each do |i|
          kmer = (seq[i..(i+k_size-1)]).upcase
          kmers << kmer
        end
        kmers.each_with_index do |kmer, index|
          @graph.add(kmer, id, index)
          @graph.add_edge(kmers[index-1], kmer) if index > 0
        end
      end
      set = Set.new
      set << 0
      set << 1
      start = @graph.first_node set

      assert_equal "TTGGAATCGGTGACCGGCATGAATTTGACAG", start

    end

    should "make one sequence out of three that overlap" do
      seq1 = "TAGCGGCCACTGAAAACTAGAATTTCCACCAAAGTTCACGAAGAGCGCGCGACTCATTCACCGCGAAGACTCACTTCGGTTTAGCGGA"
      seq2 = "TCACCGCGAAGACTCACTTCGGTTTAGCGGATGTTCACACCAATTAATGCTGCGTCCTATTGGTTTCTAGCCCATACGGCGCATACATACATACGGT"
      seq3 = "TTGGTTTCTAGCCCATACGGCGCATACATACATACGGTCCGGGATTCCATCCCACGATAGAAGGAGTCCGGAGTGCTCTATCTA"

      k_size = 31
      [seq1, seq2, seq3].each_with_index do |seq, id|
        kmers = []
        (0..seq.length-k_size).each do |i|
          kmer = (seq[i..(i+k_size-1)]).upcase
          kmers << kmer
        end
        kmers.each_with_index do |kmer, index|
          @graph.add(kmer, id, index)
          @graph.add_edge(kmers[index-1], kmer) if index > 0
        end
      end
      set = Set.new
      # set << 0
      set << 1
      set << 2
      start = @graph.first_node set

      assert_equal "TAGCGGCCACTGAAAACTAGAATTTCCACCA", start

    end

    should "check that node has contigs" do
      node = Clustemble::Node.new "ACGT"
      node.add_contig(1,1)
      node.add_contig(2,1)
      node.add_contig(3,1)
      assert_equal true, node.has_all_contigs?([1,2,3])
      assert_equal true, node.has_all_contigs?([1,2])
      assert_equal false, node.has_all_contigs?([1,2,4])
    end

    should "remove contig info" do
      node = Clustemble::Node.new "ACGT"
      node.add_contig(1,2)
      node.add_contig(2,5)
      node.add_contig(3,10)
      assert_equal true, node.has_all_contigs?([1,2,3])
      node.remove 3
      assert_equal false, node.has_all_contigs?([1,2,3])
      assert_equal true, node.has_all_contigs?([1,2])
    end

    should "rename contig on node" do
      node = Clustemble::Node.new "ACGT"
      node.add_contig(1,2)
      node.rename(1,2)
      assert_equal false, node.has_contig?(1)
      assert_equal true, node.has_contig?(2)
    end

    should "add to index" do
      node = Clustemble::Node.new "ACGT"
      node.add_contig(1,2)
      node.add_contig(2,3)
      node.add_to_index(1,2)
      assert_equal 4, node.contigs[0][:index], "index"

    end


  end

end
