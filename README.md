# Projet S1

Projet : développement d'une interface web pour la sélection positive

## Installation

### Prérequis

- [Python](https://www.python.org/downloads/) ainsi que les paquets suivants :
	- Bio
	- ete3
	- lxml
	- docopt
- [Node.js](https://nodejs.org/fr/)
	- npm

### Cloner le dépôt

`git clone git@github.com:gennaBnh/Positive-Selection.git`

## Utilisation

### Générer l'arbre à partir de l'exemple

Le script *example.sh* construira l'arbre XML correspondant aux données d'exemple fournies.

`python3 genere_xml.py <fichier_arbre> <fichier_alignement> <fichier_résultats>`

`cp <fichier_arbre>.xml input_tree.xml`

### Lancer le serveur

`npm start`

Par défaut, l'interface sera accessible à l'adresse [http://localhost:8888/](http://localhost:8888/), définie dans le fichier *selpos.js* :

```
const PORT = 8888;
const HOST = '127.0.0.1';
```
