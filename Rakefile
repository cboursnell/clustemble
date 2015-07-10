require 'rake/testtask'

Rake::TestTask.new do |t|
  t.libs << 'test'
end

desc "Run tests"
task :default => :test

Rake::TestTask.new do |t|
  t.name = :adj
  t.libs << 'test'
  t.test_files = ['test/test_adjacencylist.rb']
end

Rake::TestTask.new do |t|
  t.name = :clust
  t.libs << 'test'
  t.test_files = ['test/test_clustemble.rb']
end
