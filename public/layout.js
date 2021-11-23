d3.layout = {};

d3.layout.dendrogram = function () {
  var layout = d3.tree();
  var layoutDendrogram = function (root) {
    layout(root);
    root.each(function (d) {
      var temp = d.x;
      d.x = d.y;
      d.y = temp;
    });
  }
  layoutDendrogram.nodeSize = function (x) {
    var temp = x[0];
    x[0] = x[1];
    x[1] = temp;
    layout = layout.nodeSize(x);
    return layoutDendrogram;
  }
  layoutDendrogram.Size = function (x) {
    layout = layout.Size(x);
    return layoutDendrogram;
  }
  return layoutDendrogram;
}

d3.layout.cladogram = function () {
  var layout = d3.cluster();
  var layoutCladogram = function (root) {
    layout(root);
    root.each(function (d) {
      var temp = d.x;
      d.x = d.y;
      d.y = temp;
    });
  }
  layoutCladogram.nodeSize = function (x) {
    var temp = x[0];
    x[0] = x[1];
    x[1] = temp;
    layout = layout.nodeSize(x);
    return layoutCladogram;
  }
  layoutCladogram.Size = function (x) {
    layout = layout.Size(x);
    return layoutCladogram;
  }
  return layoutCladogram;
}

d3.layout.radial = function () {
  var layout = d3.cluster()
    .separation(function(a, b) {
      return (a.parent == b.parent ? 1 : 1) / a.depth;
    });
  var layoutRadial = function (root) {
    layout(root);
    root.each(function (d) {
      var coor = _project(d.x,d.y)
      d.x = coor[1];
      d.y = coor[0];
    });
  }
  layoutRadial.nodeSize = function (x) {
    var temp = x[0];
    x[0] = x[1];
    x[1] = temp;
    layout = layout.nodeSize(x);
    return layoutRadial;
  }
  layoutRadial.Size = function (x) {
    layout = layout.Size(x);
    return layoutRadial;
  }
  function _project(x, y) {
    var angle = (x - 90) / 180 * Math.PI,
      radius = y;
    return [radius * Math.cos(angle), radius * Math.sin(angle)];
  }
  return layoutRadial;
}
