var express = require('express');
var router = express.Router();
var multer = require('multer');
const upload = multer({ dest: 'uploads/' });
const browser = require('browser-detect');
var bodyParser = require('body-parser');
router.use(bodyParser.urlencoded({ extended: true })); // support encoded bodies
router.use(bodyParser.json()); // support json encoded bodies
if (typeof localStorage === "undefined" || localStorage === null) {
  var LocalStorage = require('node-localstorage').LocalStorage;
  localStorage = new LocalStorage('./scratch');
}
const util = require('util')

// =====================
// Code perso
// =====================

var url = require('url');

// =====================

// =====================
// Routages GET
// =====================

// Home
// ----
router.get('/', function(req, res, next) {
  res.render('index', { title: 'HOGENOM'});
});

// GET display
// -----------
router.get('/display', function(req, res, next) {
  console.log('Accès à /display');
  const fs = require('fs');
  fs.readFile('input_tree.xml', 'utf8' , (err, data) => {
    if (err) {
      res.render('error.ejs', {message:"Erreur de lecture",error:err});
    }
    var xml_digester = require("xml-digester");
    var handler = new xml_digester.OrderedElementsHandler("eventType");
    var options = {
      "handler": [{
        "path": "eventsRec/*",
        "handler": handler
      }]
    };
    var digester = xml_digester.XmlDigester(options);
    digester.digest(data, function(err, results) {
      if (err) {
        console.log(err);
        return
      }
      var JSONtree = JSON.stringify(results);
      // Séquence à mettre en valeur :
      var JSONpattern = JSON.stringify("0:NM_001193307dot1_hom_Sap_SAMD9");
      res.render('displaytree.ejs', {arbre:JSONtree,pattern:JSONpattern});
    });
  });
});

// GET taxodico
// ------------
router.get('/taxodico', function(req, res, next) {
  console.log('Accès à /taxodico');
  res.render('gettaxojson.ejs', {species:JSON.stringify([]),colour:JSON.stringify([])});
});

// =====================
// Routages POST
// =====================
router.post('/display', upload.single('file'), function(req, res, next) {
  console.log('Accès à /display');
  const fs = require('fs');
  const fname = 'uploads/' + req.file.filename
  fs.readFile(fname, 'utf8' , (err, data) => {
    if (err) {
      res.render('error.ejs', {message:"Erreur de lecture",error:err});
    }
    var xml_digester = require("xml-digester");
    var handler = new xml_digester.OrderedElementsHandler("eventType");
    var options = {
      "handler": [{
        "path": "eventsRec/*",
        "handler": handler
      }]
    };
    var digester = xml_digester.XmlDigester(options);
    digester.digest(data, function(err, results) {
      if (err) {
        console.log(err);
        return
      }
      var JSONtree = JSON.stringify(results);
      // Séquence à mettre en valeur :
      var JSONpattern = JSON.stringify("0:NM_001193307dot1_hom_Sap_SAMD9");
      res.render('displaytree.ejs', {arbre:JSONtree,pattern:JSONpattern});
    });
  });
});


// ================
// Début modifs
// ================

// GET test-affichage
// -----------
router.get('/test-affichage-graphe', function(req, res, next) {
  var pathname = url.parse(req.url).pathname;
  console.log('Accès à ' + pathname);
  res.render('test_affichage_graphe.ejs');
});


// GET test-formulaire
// -----------
router.get('/test-formulaire', function(req, res, next) {
  var pathname = url.parse(req.url).pathname;
  console.log('Accès à ' + pathname);
  res.render('test_formulaire.ejs');
});

// GET test-affichage-fichier
// -----------
router.get('/test-affichage-fichier', function(req, res, next) {
  var pathname = url.parse(req.url).pathname;
  console.log('Accès à ' + pathname);

  fs.readFile('input_tree.xml', 'utf8' , (err, data) => {
    if (err) {
      res.render('error.ejs', {message:"Erreur de lecture",error:err});
    }
    var xml_digester = require("xml-digester");
    var handler = new xml_digester.OrderedElementsHandler("eventType");
    var options = {
      "handler": [{
        "path": "eventsRec/*",
        "handler": handler
      }]
    };
    var digester = xml_digester.XmlDigester(options);
    digester.digest(data, function(err, results) {
      if (err) {
        console.log(err);
        return
      }
      var JSONtree = JSON.stringify(results);
      // Séquence à mettre en valeur :
      var JSONpattern = JSON.stringify("0:NM_001193307dot1_hom_Sap_SAMD9");
      res.render('displaytree.ejs', {arbre:JSONtree,pattern:JSONpattern});
    });
  });

  // var mydata = ''
  // const fs = require('fs');
  // fs.readFile('routes/uploads/upload-fichier-test-formulaire.txt', 'utf8', (err, data) => {
  //   if (err) {
  //     console.log('Error', err);
  //   }
  //   if (data) {
  //     mydata = data;
  //     var stat_results = JSON.stringify(mydata);
  //     console.log('--------- Results: ' + stat_results);
  //     res.render('test_affichage_fichier.ejs', {resultats:stat_results});
  //   }
  //   console.log('Data: ' + data);
  //   console.log('Mydata: ' + mydata);
  // })
  // // var mavar = JSON.parse('[[1,0.2],[2,0.33],[3,1.0]]');
  // // var stat_results = JSON.stringify('1:0.2 2:0.33 3:1.0 4:0.33 5:0.2');
});

const formidable = require('formidable');
// POST submit-test
// -----------
router.post('/submit-test', function(req, res) {
  var pathname = url.parse(req.url).pathname;
  console.log('Accès à ' + pathname);
  // console.log(req);
  // console.log('texte : '+req.body.u_texte);
  new formidable.IncomingForm().parse(req)
    .on('field', (name, field) => {
      console.log('Field', name, field);
    })
    .on('fileBegin', (name, file) => {
      console.log('Old path:', file.path);
      file.name = 'upload-' + file.name;
      file.path = __dirname + '/uploads/' + file.name;
      console.log('New name:', file.name);
      console.log('New path:', file.path);
    })
    .on('file', (name, file) => {
      // console.log('File', name, file);
      console.log('File', file.name, 'received, here is some information about it:');
      console.log('Size:', file.size, 'bytes');
      console.log('Type:', file.type);
      console.log('Path:', file.path)
    })
    .on('aborted', () => {
      console.log('File upload aborted by the user');
    })
    .on('error', (err) => {
      console.log('Error', err);
      throw err;
    })
    .on('end', () => {
      res.end();
    });
});

// ================
// Fin modifs
// ================

module.exports = router;
