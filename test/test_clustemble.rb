#!/usr/bin/env	ruby

require 'helper'

class TestClustemble < Test::Unit::TestCase

  context 'clustemble' do

    setup do
      @clust = Clustemble::Clustemble.new 31
    end

    teardown do
    end

    should "do something" do
      file = File.join(File.dirname(__FILE__), 'data', 'cluster.fa')
      @clust.add_fasta file
    end

    should "kmerise sequence" do
      list = @clust.kmerise "TTGGAATCGGTGACCGGCATGAATTTGACAGAACTCGAGGCGATT"
      assert_equal 15, list.size
      assert_equal 31, list[0].length
    end

    should "find start" do
      file = File.join(File.dirname(__FILE__), 'data', 'cluster.fa')
      @clust.add_fasta file
      assert_equal "CACAAAACTAAGATCTTGTTCATTTCCTATG", @clust.find_start
    end

    should "traverse" do
      file = File.join(File.dirname(__FILE__), 'data', 'cluster.fa')
      @clust.add_fasta file
      start = @clust.find_start
      seq = @clust.traverse start
      File.open("cluster.aln.fa", "wb") do |out|
        seq.each_with_index do |str,i|
          out.write ">contig#{i}\n"
          out.write "#{str[0..-2]}\n"
        end
      end
      # assert_equal 3, seq.length
    end

    should "traverse a simple graph" do
      clust = Clustemble::Clustemble.new 15
      file = File.join(File.dirname(__FILE__), 'data', 'test.fa')
      clust.add_fasta file
      start = clust.find_start
      seq = clust.traverse start
      seq.each_with_index do |str,i|
        puts ">contig#{i}"
        puts str[0..-2]
      end
      assert_equal "ACATACAGGCAATAAGCAACACC", seq[0][-24..-2]
    end

    should "traverse a slightly more complex graph" do
      clust = Clustemble::Clustemble.new 15
      file = File.join(File.dirname(__FILE__), 'data', 'test2.fa')
      clust.add_fasta file
      start = clust.find_start
      seq = clust.traverse start
      seq.each_with_index do |str,i|
        puts ">contig#{i}"
        puts str[0..-2]
      end
      assert_equal 1, seq.length
      assert_equal 968, seq[0].length
    end

    should "traverse an even more complex graph" do
      clust = Clustemble::Clustemble.new 15
      file = File.join(File.dirname(__FILE__), 'data', 'test3.fa')
      clust.add_fasta file
      start = clust.find_start
      seq = clust.traverse start
      seq.each_with_index do |str,i|
        puts ">contig#{i}"
        puts str[0..-2]
      end
      assert_equal 2, seq.length
      assert_equal 968, seq[0].length
      assert_equal 618, seq[1].length
    end

    should "traverse a graph with two starts" do
      clust = Clustemble::Clustemble.new 15
      file = File.join(File.dirname(__FILE__), 'data', 'test4.fa')
      clust.add_fasta file
      start = clust.find_start
      seq = clust.traverse start
      seq.each_with_index do |str,i|
        puts ">contig#{i}"
        puts str[0..-2]
      end
      assert_equal 2, seq.length
      assert_equal 968, seq[0].length
      assert_equal 618, seq[1].length
    end

  end

end
