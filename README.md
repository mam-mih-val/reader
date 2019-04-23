# reader

This project is made for flow analysis in ultrarelativistic heavy ions collisions. The inteface is build for DataTree format only.

1. Build

Firstly, 

    $ git clone  https://github.com/mam-mih-val/reader

Then,

    $ cd reader
    $ mkdir build

You need the last version of DataTree (hades branch), which you can get there, follow the instruction in README https://gitlab.cern.ch/na61-hic/DataTree

After the DataTree installation, make

    $ export DATATREE_HOME=/path/to/data/tree/installation
    $ cmake ..
    $ make
    
2. Using

Change current directory to your build and than simply launch DT_Reader as bash-script:

    $ ./DT_Reader
