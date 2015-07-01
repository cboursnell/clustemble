# coding: utf-8
require File.expand_path('../lib/clustemble/version', __FILE__)

Gem::Specification.new do |gem|
  gem.name          = "clustemble"
  gem.version       = Clustemble::VERSION::STRING.dup
  gem.authors       = ["Chris Boursnell"]
  gem.email         = ["cmb211@cam.ac.uk"]
  gem.summary       = %q{TODO: Write a short summary. Required.}
  gem.description   = %q{TODO: Write a longer description. Optional.}
  gem.homepage      = ""
  gem.license       = "MIT"

  gem.files       = `git ls-files`.split("\n")
  gem.executables = ["clustemble"]
  gem.require_paths = ["lib"]

  gem.add_dependency 'trollop','~> 2.1', '>= 2.1.1'
  gem.add_dependency 'bio', '~> 1.4', '>= 1.4.3'
  gem.add_dependency 'rgl', '~> 0.5', '>= 0.5.0'

  gem.add_development_dependency 'rake', '~> 10.3', '>= 10.3.2'
  gem.add_development_dependency 'turn', '~> 0.9', '>= 0.9.7'
  gem.add_development_dependency 'simplecov', '~> 0.8', '>= 0.8.2'
  gem.add_development_dependency 'shoulda-context', '~> 1.2', '>= 1.2.1'
  gem.add_development_dependency 'coveralls', '~> 0.7'
end
