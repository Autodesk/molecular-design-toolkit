#!/usr/bin/env bash

rvm install 2.4.1
rvm use 2.4.1

gem install bundler travis
travis

git clone https://github.com/travis-ci/travis-build.git
cd travis-build
bundler install