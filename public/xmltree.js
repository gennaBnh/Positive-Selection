xmlparser = {};

xmlparser.parse = function(xmlstr) {
  console.log("XML=");
  console.log(xmlstr);
}

xmlparser.flatTree = function(treeRoot,config = {}) {
  var virtualRoot = {
    name : "Out",
    eventsRec : treeRoot.eventsRec,
    clade : [treeRoot]
  }
  treeRootNode = d3.hierarchy(virtualRoot,function(d) {
    return d.clade;
  });
  treeRootNode.each(function (node) {
    if(node.children && node.children.length) {
      node.children.forEach(function (child,posChild) {
        var listEvents = child.data.eventsRec;
        var newEvent = null;
        var startNode = undefined;
        var currentNode = undefined;
        while (listEvents && (newEvent = listEvents.shift())){
          switch (newEvent.eventType) {
            case "speciationLoss":
              var newChildName = child.data.name+"_SpL";
              var lossChildName = child.data.name+"_Loss";
              if(config.speciationLoss == false) {
                var newChild = createNewSubTree(newChildName,newEvent);
              }
              else {
                var newChild = createNewSubTreeWithChild(newChildName,newEvent,lossChildName,"undefined");
              }
              if(!startNode && !currentNode) {
                startNode = newChild;
                currentNode = newChild;
              } else {
                currentNode.clade.push(newChild);
                currentNode = newChild;
              }
              break;
            case "speciationOutLoss":
              var newChildName = child.data.name+"_SpOL";
              var lossChildName = child.data.name+"_Loss";
              if(config.speciationOutLoss == false) {
                var newChild = createNewSubTree(newChildName,newEvent);
              }
              else {
                var newChild = createNewSubTreeWithChild(newChildName,newEvent,lossChildName,newEvent.speciesLocation);
              }
              if(!startNode && !currentNode) {
                startNode = newChild;
                currentNode = newChild;
              }
              else {
                currentNode.clade.push(newChild);
                currentNode = newChild;
              }
              break;
            case "transferBack":
              var newChildName = child.data.name+"_TrB";
              var lossChildName = child.data.name+"_Loss";
              if(config.transferBack == false) {
                var newChild = createNewSubTree(newChildName,newEvent);
              }
              else {
                var newChild = createNewSubTreeWithChild(newChildName,newEvent,lossChildName,"out");
              }
              if(!startNode && !currentNode) {
                startNode = newChild;
                currentNode = newChild;
              }
              else {
                currentNode.clade.push(newChild);
                currentNode = newChild;
              }
              break;
            default:
              child.data.eventsRec = [newEvent];
          }
        }
        if(startNode && currentNode) {
          node.data.clade[posChild] = startNode;
          currentNode.clade.push(child.data);
        }
      });
    }
  });
  return virtualRoot.clade[0];
}

function createNewSubTree(nodeName,nodeEvent) {
  return {
    name : nodeName,
    eventsRec : [nodeEvent],
    clade : []
  }
}

function createNewSubTreeWithChild(nodeName,nodeEvent,childName,childSpeciesLocation) {
  return {
    name : nodeName,
    eventsRec : [nodeEvent],
    clade : [
      {
        name: childName,
        eventsRec : [{eventType: 'loss' , speciesLocation: childSpeciesLocation}]
      }
    ]
  }
}
