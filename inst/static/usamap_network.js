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

// filter out NA values from location (done in the main function)
// added here so it has scope outside of the onclikc function
var pointdatafiltered = [];


makelegend = function(maxval){
    // first f=r=emvoe old legends
    d3.selectAll(".legend").remove();

    var radius = d3.scale.sqrt()
    .domain([0, maxval])
    .range([0, 20]);

    var legend =d3.select("svg").append("g")
        .attr("class", "legend")
        .attr("transform", "translate(" + (width - 50) + "," + (height - 20) + ")")
      .selectAll("g")
        .data([maxval/10, maxval/2, maxval])
      .enter().append("g");

    legend.append("circle")
        .attr("cy", function(d) { return -radius(d); })
        .attr("r", radius);

    legend.append("text")
        .attr("y", function(d) { return -2 * radius(d); })
        .attr("dy", "1.3em")
        .text(d3.format(".1s"));
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


drawLines = function(){
    //linedata and pointdata are included in the global scope
    links = []

    for(var i=0, len=linedata.length-1; i<len; i++){
    // (note: loop until length - 1 since we're getting the next
    //  item with i+1)
        var source = linedata[i].source;
        var target = linedata[i].target;
        var distancevalue = linedata[i].value;
        //console.log(source);
        //console.log(target);
        //console.log(distancevalue);

        loc1 = pointdata.find(function(d){return d.samplename === source})
        loc2 = pointdata.find(function(d){return d.samplename === target})

        links.push({
            type: "LineString",
            coordinates: [
                [ loc1.Longitude, loc1.Latitude ],
                [ loc2.Longitude, loc2.Latitude ]
            ]
        });

    var pathArcs = d3.select("svg")
        .append("g")
        .selectAll(".arc")
        .data(links);

    //enter
    pathArcs.enter()
        .append("path").attr({
            'class': 'arc'
        }).style({
            fill: 'none',
        });

    //update
    pathArcs.attr({
            //d is the points attribute for this path, we'll draw
            //  an arc between the points using the arc function
            d: path
        })
        .style({
            stroke: '#0000ff',
            'stroke-width': '2px'
        })
    }
}

drawMap = function() {
    // Load BaseMap is from a different file, allowing the basemap to change from USA to World, etc/
    load_basemap()

    //draw the lines first
    drawLines()

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


}






