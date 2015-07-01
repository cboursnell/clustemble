require 'rake/testtask'
require 'rake/extensiontask'
require 'bundler/setup'

Rake::ExtensionTask.new('transrate') do |ext|
  ext.lib_dir = "lib/transrate"
end

Rake::TestTask.new do |t|
  t.libs << 'test'
end
