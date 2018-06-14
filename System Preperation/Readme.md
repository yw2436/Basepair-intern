The test server I used is a Microsoft Azure server(CentOS 7). I need to do a lot preparation before using it.
During installing all the pipelines, I installed some system packages that are needed using 'sudo yum install'

List:
  gcc-c++
  glibc-static
  ncurses-devel
  zlib-devel
  bzip-devel
  bzip2-devel
  xz-devel
  libcurl-devel.x86_64
  python-setuptools.noarch
  git
  python-pip
And using pip, I installed some python packages: numpy, methylpy
