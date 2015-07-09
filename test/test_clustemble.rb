#!/usr/bin/env	ruby

require 'helper'

class TestClustemble < Test::Unit::TestCase

  context 'clustemble' do

    setup do
      @clust = Clustemble::Clustemble.new 31
    end

    teardown do
    end

    # should "kmerise sequence" do
    #   list = @clust.kmerise "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   assert_equal 15, list.size
    #   assert_equal 31, list[0].length
    # end

    # should "kmerise sequence with lower case characters" do
    #   list1 = @clust.kmerise "TTGGAATCGGTGACCGGcATGAATTTGACAGAACTCGAGGCGATT"
    #   list2 = @clust.kmerise "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   assert_equal 15, list1.size
    #   assert_equal 31, list1[0].length
    #   assert list1==list2
    # end

    # should "add sequence" do
    #   @clust.add_seq 0, "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   assert_equal 15, @clust.graph.size
    #   assert_equal 1, @clust.graph.edges["TTGGAATCGGTGACCGGCATGAATTTGACAG"].size
    #   assert_equal 14, @clust.graph.num_edges
    # end

    # should "add fasta file" do
    #   file = File.join(File.dirname(__FILE__), 'data', 'test.fa')
    #   @clust.add_fasta file
    #   assert_equal 937, @clust.graph.size
    # end

    # should "find sequence is redundant" do
    #   seq1 ="TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   seq2 = seq1[4..41]
    #   @clust.add_seq 0, seq1
    #   @clust.add_seq 1, seq2
    #   assert_equal 15, @clust.graph.size
    #   @clust.extract_seqs
    # end

    should "combine overlapping sequences" do
      seq0 = "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      seq1 ="TCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      seq2 ="TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAG"
      @clust.add_seq 0, seq1
      @clust.add_seq 1, seq2
      assert_equal 15, @clust.graph.size
      seqs = @clust.extract_seqs
      assert_equal 0, seqs.keys.first
      assert_equal seq0, seqs.values.first
    end

    should "combine three overlapping sequences " do
      seq3 = "CACAAAACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTATAAGCAAATGGCAAA"
      seq3 << "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
      seq1 = "CACAAAACTAAGATCTTGTTCATTTCCTATGAC"
      seq2 = "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
      @clust.add_seq 1, seq1
      assert_equal 3, @clust.graph.size
      @clust.add_seq 2, seq2
      assert_equal 8, @clust.graph.size
      assert_equal 2, @clust.graph.starts.size
      @clust.add_seq 3, seq3
      assert_equal 1, @clust.graph.starts.size
      assert_equal 66, @clust.graph.size
      p @clust.extract_seqs

    end

    # should "find starts" do
    #   file = File.join(File.dirname(__FILE__), 'data', 'cluster.fa')
    #   @clust.add_fasta file
    #   list = @clust.find_starts
    #   assert_equal "CACAAAACTAAGATCTTGTTCATTTCCTATG", list[0]
    # end


    # should "keep sequence because it's a different isoform" do
    #   @clust.add_seq 0, "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   @clust.add_seq 1, "TTGGAATCGGTGACCGGCATTTTGACAGAACTCGAGGCGATT"
    #   assert_equal 27, @clust.graph.size, "graph size"
    #   @clust.align_seq(1, "TTGGAATCGGTGACCGGCATTTTGACAGAACTCGAGGCGATT")
    # end

  end

end
