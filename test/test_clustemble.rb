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
    #   assert_equal 15, list.size, "list length"
    #   assert_equal 31, list[0].length, "kmer length"
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
    #   # puts "ADDING 0: #{seq1}"
    #   @clust.add_seq 0, seq1
    #   # puts "ADDING 1: #{seq2}"
    #   @clust.add_seq 1, seq2
    #   assert_equal 17, @clust.graph.size
    #   assert_equal seq1, @clust.extract_seqs[100]
    # end

    # should "deal with the same kmer appearing twice in a sequence" do
    #   seq1 = "TTGGAATCGGTGACCGGCATGAATTTGACAGATT"
    #   seq2 = "AATTGGAATCGGTGACCGGCATGAATTTGACAGATTGGAATCGGTGACCGGCATGAATTTGACAGTA"
    #   list = @clust.kmerise 0, seq2
    #   assert_equal 38, list.size, "list size"
    #   hash={}
    #   list.each do |i|
    #     hash[i]||=0
    #     hash[i] += 1
    #   end
    #   assert_equal 36, hash.size, "hash size"
    #   @clust.add_seq 0, seq1
    #   @clust.add_seq 1, seq2
    #   seqs = @clust.extract_seqs
    #   assert seqs.size > 0, "seq hash size"
    #   assert_equal seq2, seqs.to_a.first.last, "seq2"
    # end

    # should "deal with some kmers appearing multiple times" do
    #   seq1 = "GGCGCCATACCCATCTTACGTAATGTCAATAAAAACATGG"
    #   seq2 = "GGCGCCATACCCATCTTACGTAATGTCAATAAAAACATGGTTTCCATACCCATCTTACGTAA"
    #   seq2 << "TGTCAATAAAAACATGGCCCCCATACCCATCTTACGTAATGTCAATAAAAACATG"
    #   @clust.add_seq 1, seq1
    #   @clust.add_seq 2, seq2
    #   seqs = @clust.extract_seqs
    #   assert_equal seq2, seqs.to_a.first.last, "seq2"
    # end

    # should "deal with some kmers appearing multiple times again" do
    #   seq1 = "AATTTAGTGGTCTAGAACACCTTATTAGTGCCAGCGGTTTAACCATTGACGGTAGGAGTTGG"
    #   seq1 << "ACGTTCCATGCTGGTAGATCTGCCAGCGGTTTAACCATTGACGGTAGGAGTAG"
    #   seq2 = "TGCCAGCGGTTTAACCATTGACGGTAGGAGTAGCTCCCGACACCTCTGAGCAGCTTACGACGA"
    #   seq2 << "CTAATAGATAGGCCTCGGCGCTTATCAGGAATGTGGCGTCCAAAGTGAATGTGCCAGCGGTT"
    #   seq2 << "TAACCATTGACGGTAGGAGTCGCGTTCACGATGGGTCACCCTAATAGGATTACCA"
    #   res = "AATTTAGTGGTCTAGAACACCTTATTAGTGCCAGCGGTTTAACCATTGACGGTAGGAGTTGGAC"
    #   res << "GTTCCATGCTGGTAGATCTGCCAGCGGTTTAACCATTGACGGTAGGAGTAGCTCCCGACACCT"
    #   res << "CTGAGCAGCTTACGACGACTAATAGATAGGCCTCGGCGCTTATCAGGAATGTGGCGTCCAAAG"
    #   res << "TGAATGTGCCAGCGGTTTAACCATTGACGGTAGGAGTCGCGTTCACGATGGGTCACCCTAATA"
    #   res << "GGATTACCA"
    #   @clust.add_seq 1, seq1
    #   @clust.add_seq 2, seq2
    #   seqs = @clust.extract_seqs
    #   assert_equal res, seqs.to_a.first.last, "seq"
    # end

    # should "make one sequence out of three that overlap" do
    #   seq1 = "TAGCGGCCACTGAAAACTAGAATTTCCACCAAAGTTCACGAAGAGCGCGCGACTCATTCACCGCGAAGACTCACTTCGGTTTAGCGGA"
    #   seq2 = "TCACCGCGAAGACTCACTTCGGTTTAGCGGATGTTCACACCAATTAATGCTGCGTCCTATTGGTTTCTAGCCCATACGGCGCATACATACATACGGT"
    #   seq3 = "TTGGTTTCTAGCCCATACGGCGCATACATACATACGGTCCGGGATTCCATCCCACGATAGAAGGAGTCCGGAGTGCTCTATCTA"
    #   @clust.add_seq 1, seq1
    #   @clust.add_seq 2, seq2
    #   @clust.add_seq 3, seq3
    #   seqs = @clust.extract_seqs
    #   ans = "TAGCGGCCACTGAAAACTAGAATTTCCACCAAAGTTCACGAAGAGCGCGCGACTCATTCACCGC"
    #   ans << "GAAGACTCACTTCGGTTTAGCGGATGTTCACACCAATTAATGCTGCGTCCTATTGGTTTCTAG"
    #   ans << "CCCATACGGCGCATACATACATACGGTCCGGGATTCCATCCCACGATAGAAGGAGTCCGGAGT"
    #   ans << "GCTCTATCTA"
    #   assert_equal ans, seqs.to_a.first.last, "seq"
    # end

    # should "combine overlapping sequences" do
    #   # ********----
    #   # ----********
    #   #             v-----------------seq1----------------v
    #   seq0 = "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   #       ^--------------seq2-------------------^
    #   seq1 ="TCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
    #   seq2 ="TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAG"
    #   @clust.add_seq 1, seq1
    #   @clust.add_seq 2, seq2
    #   assert_equal 17, @clust.graph.size
    #   seqs = @clust.extract_seqs
    #   assert_equal seq0, seqs.to_a.first.last, "sequence"
    #   # assert_equal 1, seqs.keys.first
    #   # assert_equal seq0, seqs.values.first
    # end

    # should "combine three overlapping sequences " do
    #   seq3 = "CACAAAACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTATAAGCAAATGGCAAA"
    #   seq3 << "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
    #   seq1 = "CACAAAACTAAGATCTTGTTCATTTCCTATGAC"
    #   seq2 = "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
    #   @clust.add_seq 1, seq1
    #   assert_equal 4, @clust.graph.size, "graph size"
    #   @clust.add_seq 2, seq2
    #   p @clust.graph
    #   assert_equal 10, @clust.graph.size, "graph size"
    #   assert_equal 2, @clust.graph.starts.size, "graph start size"
    #   @clust.add_seq 3, seq3
    #   assert_equal 1, @clust.graph.starts.size, "graph start size"
    #   assert_equal 69, @clust.graph.size, "graph size"
    #   result = @clust.extract_seqs
    #   assert_equal seq3, result.to_a.first.last
    # end

    should "return two sequences" do
      seq1 =  "CACAAAACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTA"
      seq1 << "TAAGCAAATGGCAAA"
      seq1 << "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
      seq2 =  "CACAAAACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTA"
      seq2 << "CTATATTTATCATCAAATGTCAAGATGCATCTAAT"
      seq3 =       "AACTAAGATCTTGTTCATTTCCTATGACAATAACATTATTA"
      seq3 << "TAAGCA"

      # 3 ----***************************-----------------------------
      # 1 ************************************************************
      # 2 ****************************----------**********************

      # puts "TEST: adding contig id 1"
      @clust.add_seq 0, seq1
      assert_equal 67, @clust.graph.size, "graph size after adding seq1"
      # puts "TEST: adding contig id 2"
      @clust.add_seq 1, seq2
      assert_equal 97, @clust.graph.size, "graph size after adding seq2"
      assert_equal 2, @clust.extract_seqs.size, "extracted seq size"

      # puts "TEST: adding contig id 3"
      @clust.add_seq 2, seq3
      seqs = @clust.extract_seqs
      p seqs
      assert_equal 2, seqs.size
      assert_equal seq1, seqs.to_a[0].last, "seq1"
      assert_equal seq2, seqs.to_a[1].last, "seq2"
    end

    # should "run on an actual fasta file" do
    #   file = File.join(File.dirname(__FILE__), 'data', 'cluster.fa')
    #   @clust.add_fasta file
    #   seqs = @clust.extract_seqs
    #   assert_equal 4, seqs.size
    # end

    # should "run on another actual fasta file" do
    #   file = File.join(File.dirname(__FILE__), 'data', 'test4.fa')
    #   @clust.add_fasta file
    #   seqs = @clust.extract_seqs
    #   assert_equal 2, seqs.size
    # end

    # should "return consensus fasta file" do
    #   file = File.join(File.dirname(__FILE__), 'data', 'test4.fa')
    #   fasta = @clust.consensus file, "test"
    #   assert_equal 1385, fasta.length
    # end

  end

end
