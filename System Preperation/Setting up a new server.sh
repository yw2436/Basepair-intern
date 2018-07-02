install pip:
  sudo yum -y update                          #not root usr, use sudo.
  sudo yum --enablerepo=extras install epel-release
  sudo yum search python2-pip                 #the name may change, it was python-pip, and now it's python2-pip. May change again.
  sudo yum install #package name#
  sudo pip install --upgrade pip              #upgrade
install methylpy(automaticlly including numpy, scipy):
  sudo pip install methylpy
install sratools:
  #first copy the binary files downloaded from their website to the server
  #and move them to the bin directory, and add them to PATH. Also, change the permissions of the programs under the bin directory
make the /mnt/resource directory usable:
  sudo chown yuhan /mnt/resource
install bowtie/bowtie2:


  
 
