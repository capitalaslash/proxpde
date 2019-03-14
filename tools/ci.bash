#! /bin/bash

parallel=12
ci_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
envfile=$ci_dir/ci.env
logfile=$ci_dir/ci.log
lockfile=$ci_dir/.ci.lock

tmp_dir=/tmp/ci

message () {
  echo -e "\e[31m$1\e[0m"
}

quit () {
  echo $1 | tee -a $logfile
  rm $lockfile
  exit $2
}

if [ -f $lockfile ]; then
  echo "ci in use, wait for it to end or remove lock file $lockfile"
  exit 1
fi
touch $lockfile

echo "execution start: $(date)" | tee $logfile
echo "commit: $(git rev-parse HEAD)" | tee -a $logfile
mkdir -p $tmp_dir
pushd $tmp_dir || quit "tmp dir $tmp_dir not found" 2

message "cloning..."
git clone /home/cervone/repo/minifem >> $logfile || quit "clone failed in $tmp_dir" 3
pushd minifem/

message "loading environment..."
source $envfile
module list 2>> $logfile

for build in Release Debug; do
  mkdir build
  pushd build/

  message "configuring $build build..."
  cmake .. -DCMAKE_BUILD_TYPE=$build >> $logfile || quit "cmake failed" 4

  message "compiling..."
  make -j$parallel >> $logfile || quit "make failed" 5

  pushd test/
  message "testing..."
  ctest -j$parallel >> $logfile || quit "test failed" 6
  popd

  popd
  rm -r build || quit "build rm failed" 7

  echo | tee -a $logfile
  echo | tee -a $logfile
done

popd
message "cleaning up..."
rm -r minifem || quit "clean up failed" 8
popd

quit "execution end: $(date)" 0

