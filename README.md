# Projet S1

Project : development of an interface for positive selection

## Installation

### Prerequisites

- [Python 3](https://www.python.org/downloads/)
	- Bio
	- ete3
	- lxml
- [Node.js](https://nodejs.org/fr/)
	- npm

### Cloning the repository

`git clone git@github.com:gennaBnh/Positive-Selection.git`

## Usage

### Generating an XML file for the interface

In order to create the XML file the interface needs, use the command

`python3 genere_xml.py -t <tree_file> -a <alignment_file> -s <results_file> -c <results_column>`

- You may choose which column of results should be considered when creating the XML file.
It will be set to 1 by default (0-based) in order to skip a line header, but you can change this
with the `-c`/`--col` option.
- In case some positions have no result associated with them, a replacement value of -1 will be
applied by default. This value can be changed with the `-n`/`--nostat` option.
- The output XML can be renamed with the `-o`/`--output` option.

#### Generating an XML file for the interface - example : 

An example with the files in the folder 'data' : 

`python3 genere_xml.py -t data/example_tree.tree -a data/example_alignment.fasta -s data/example_results.info -c 8`

### Launching the server

Simply use `npm start` and the interface will be available at the configured address.

The address ([http://127.0.0.1:8888/](http://localhost:8888/) by default) is defined in *selpos.js*.

```
const PORT = 8888;
const HOST = '127.0.0.1';
```
