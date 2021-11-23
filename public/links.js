var svgLinks = {};

svgLinks.shoulder = function(d) {
  var path = d3.path();
  path.moveTo(d.source.x, d.source.y);
  path.lineTo(d.source.x, d.target.y);
  path.lineTo(d.target.x, d.target.y);
  return path.toString();
}

svgLinks.line = function(d) {
  var path = d3.path();
  path.moveTo(d.source.x, d.source.y);
  path.lineTo(d.target.x, d.target.y);
  return path.toString();
}

svgLinks.radial = function(d) {
  var path = "";
  path += "M" + [d.source.x, d.source.y];

  path += "C" + [d.target.x, (d.source.y + d.target.y) / 2] + " "
  path += [d.target.x, d.target.y]+ " "
  path += [d.target.x, d.target.y];
  return path;
};
