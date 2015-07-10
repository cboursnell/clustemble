#!/usr/bin/env	ruby

require 'helper'

class TestClustemble < Test::Unit::TestCase

  context 'clustemble' do

    setup do
      @clust = Clustemble::Clustemble.new 31
    end

    teardown do
    end

    should "kmerise sequence" do
      list = @clust.kmerise "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      assert_equal 15, list.size
      assert_equal 31, list[0].length
    end

    should "kmerise sequence with lower case characters" do
      list1 = @clust.kmerise "TTGGAATCGGTGACCGGcATGAATTTGACAGAACTCGAGGCGATT"
      list2 = @clust.kmerise "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      assert_equal 15, list1.size
      assert_equal 31, list1[0].length
      assert list1==list2
    end

    should "add sequence" do
      @clust.add_seq 0, "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      assert_equal 15, @clust.graph.size
      assert_equal 1, @clust.graph.edges["TTGGAATCGGTGACCGGCATGAATTTGACAG"].size
      assert_equal 14, @clust.graph.num_edges
    end

    should "add fasta file" do
      file = File.join(File.dirname(__FILE__), 'data', 'test.fa')
      @clust.add_fasta file
      assert_equal 937, @clust.graph.size
    end

    should "find sequence is redundant" do
      seq1 ="TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      seq2 = seq1[4..41]
      @clust.add_seq 0, seq1
      @clust.add_seq 1, seq2
      assert_equal 15, @clust.graph.size
      @clust.extract_seqs
    end

    should "deal with the same kmer appearing twice in a sequence" do
      list = @clust.kmerise "TTGGAATCGGTGACCGGCATGAATTTGACAGATTGGAATCGGTGACCGGCATGAATTTGACAG"
      assert_equal 33, list.size
    end

    should "combine overlapping sequences" do
      # ********----
      # ----********
      seq0 = "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      seq1 ="TCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      seq2 ="TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAG"
      @clust.add_seq 1, seq1
      @clust.add_seq 2, seq2
      assert_equal 15, @clust.graph.size
      seqs = @clust.extract_seqs
      p seqs
      assert_equal 1, seqs.keys.first
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
      assert_equal seq3, @clust.extract_seqs[1]
    end

    should "return two sequences" do
      seq1 =  "CACAAAACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTA"
      seq1 << "TAAGCAAATGGCAAA"
      seq1 << "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
      seq2 =  "CACAAAACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTA"
      seq2 << "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
      seq3 =       "AACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTA"
      seq3 << "TAAGCA"

      # puts "TEST: adding contig id 1"
      @clust.add_seq 1, seq1
      assert_equal 66, @clust.graph.size, "graph size after adding seq1"
      # puts "TEST: adding contig id 2"
      @clust.add_seq 2, seq2
      assert_equal 95, @clust.graph.size, "graph size after adding seq2"
      assert_equal 2, @clust.extract_seqs.size, "extracted seq size"

      # puts "TEST: adding contig id 3"
      @clust.add_seq 3, seq3
      seqs = @clust.extract_seqs
      assert_equal 2, seqs.size
      assert_equal seq1, seqs[1], "seq1"
      assert_equal seq2, seqs[2], "seq2"
    end

    should "run on an actual fasta file" do
      file = File.join(File.dirname(__FILE__), 'data', 'cluster.fa')
      @clust.add_fasta file
      seqs = @clust.extract_seqs
      assert_equal 4, seqs.size
    end

    should "run on another actual fasta file" do
      file = File.join(File.dirname(__FILE__), 'data', 'test4.fa')
      @clust.add_fasta file
      seqs = @clust.extract_seqs
      assert_equal 2, seqs.size
    end

  end

end
