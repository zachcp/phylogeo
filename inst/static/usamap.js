// basic Bubble Map inspired by
// http://bl.ocks.org/mbostock/9943478
//
//

var width = 960,
    height = 500;

var projection = d3.geo.albersUsa()
    .scale(1000)
    .translate([width / 2, height / 2]);

var path = d3.geo.path()
    .projection(projection);

var selectvals = Object.keys(pointdata[0].countdata);

var select = d3.select('body')
  .append('select')
  	.attr('class','select')
    .on('change', onchange)

// filter out NA values from location (done in the main function)
// added here so it has scope outside of the onclikc function
var pointdatafiltered = [];


makelegend = function(maxval){
    // first f=r=emvoe old legends
    d3.selectAll(".legend").remove();

    var radius = d3.scale.sqrt()
    .domain([0, 0.1])
    .range([0, 20]);

    var legend = d3.select("svg").append("g")
        .attr("class", "legend")
        .attr("transform", "translate(" + (width - 50) + "," + (height - 20) + ")")
      .selectAll("g")
        .data([0.01, 0.1, 0.1])
      .enter().append("g");

    legend.append("circle")
        .attr("cy", function(d) { return -radius(d); })
        .attr("r", radius);

    legend.append("text")
        .attr("y", function(d) { return -2 * radius(d); })
        .attr("dy", "1.3em")
        .text(d3.format(".1s"));
};


onchange = function() {
    value = d3.select('select').property('value')
    maxvalue = d3.max(pointdatafiltered, function(d){return d.countdata[value]})

    var radius = d3.scale.sqrt()
        .domain([0, 0.1])
        .range([0, 20]);

    d3.selectAll("circle")
        .data(pointdatafiltered
                .sort(function(a, b) { return b.countdata[value] - a.countdata[value]; }))
        .attr("r", function(d){return radius(+d.countdata[value] / +d.totalreads )})

    makelegend(maxvalue)
};

load_basemap  = function() {
    // create the base map on SVG
    var svg = d3.select("body").append("svg")
    .attr("width", width)
    .attr("height", height);

    // create USA Basemap
    svg.insert("path", ".graticule")
      .datum(topojson.feature(usa, usa.objects.land))
      .attr("class", "land")
      .attr("d", path);

    svg.insert("path", ".graticule")
      .datum(topojson.mesh(usa, usa.objects.states, function(a, b) { return a !== b; }))
      .attr("class", "state-boundary")
      .attr("d", path);

    svg.insert("path", ".graticule")
      .datum(topojson.mesh(usa, usa.objects.states, function(a, b) { return a !== b; }))
      .attr("class", "state-boundary")
      .attr("d", path);
}


drawMap = function() {
    // Load BaseMap is from a different file, allowing the basemap to change from USA to World, etc/
    load_basemap()

    //filter points that have locations outside the map boundary
    pointdatafiltered = []

    //pointdata is supplied by the app
    //add location and pop it into the filtered array
    pointdata.forEach(function(point){
	var location = [+point.Longitude, +point.Latitude];
	var projectionlocation = projection(location)
	if (projectionlocation) {
		point.location = projectionlocation
		pointdatafiltered.push(point)
	}});


    // draw the circles
    var circles = d3.select("svg").append("svg:g")
    .attr("id", "circles");

    circles.selectAll("circle")
        .data(pointdatafiltered)
      .enter().append("svg:circle")
        .attr("cx", function(d) { return d.location[0]; })
        .attr("cy", function(d) { return d.location[1]; })
        .attr("r", 5)
        .attr("class", "bubble");

    // Add options to the select box based on the values in the
    // countdata array
   d3.select('select')
  	.attr('class','select')
    .on('change',onchange)
        .selectAll('option')
    	    .data(selectvals)
    	    .enter()
    	    .append('option')
    		.text(function (d) { return d; })
    		.attr("value", function(d) {return d});

  // add a legend
  makelegend();
}






