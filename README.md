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
    
2. Usage

Change current directory to your build and than simply launch DT_Reader as bash-script:

    $ ./DT_Reader ["input filename"] [command] [configuration] ["output filename"]

commands:

      qa -        build Quality Assurance histograms
      qvector -   build Q-vector histograms
      flow -      build Flow histograms

configiration (optional for flow & qvector)

      [1/0][1/0][1/0]
      first -     0/1 - channel selection off/on
      second -    0/1 - Using Fw-Adc, Fw-Z
      third -     0/1 - Using All spectators or only protons 
      
for example
    
    $ ./DT_Reader /path/to/input/file.root flow 001 /path/to/output/file.root
