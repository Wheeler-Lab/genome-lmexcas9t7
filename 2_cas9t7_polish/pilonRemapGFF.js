//Load gff file of reference used for polishing
var fs=require("fs");
var path=require("path");
var gff=fs.readFileSync(path.resolve(process.cwd(), process.argv[2])).toString();
gff=gff.split("\n");

//Setup data objects
var coordOffs={};
var coordMaps={};

//Setup readline for stdin input of changes pilon changes file
var readline=require("readline");
var rl=readline.createInterface({
	input: process.stdin,
	output: process.stdout,
	terminal: false
});

//On new stdin line
rl.on("line", function(line){
	var data=line.split(" ");
	//Parse cols 1 and 2 for contig ID and coordinate ranges
	for (var i=0; i<=1; i++) {
		data[i]=data[i].split(":");
		if (data[i][1].indexOf("-")!=-1) {
			data[i][1]=data[i][1].split("-");
		} else {
			data[i][1]=[data[i][1], data[i][1]];
		}
		for (var j=0; j<=1; j++) {
			data[i][1][j]=parseInt(data[i][1][j]);
		}
	}
	//Correct coordinates if change from empty
	for (var i=0; i<=1; i++) {
		if (data[i+2]!=".") data[i][1][1]++;
	}
	//Determine length change
	var dlength=(data[1][1][1]-data[1][1][0])-(data[0][1][1]-data[0][1][0]);
	if (dlength!=0) {
		//Record changes which alter length => alter coordinates
		//console.error(data[0][0], data[0][1], data[1][1]);
		//Update coordinate offsets
		if (!coordOffs[data[0][0]]) coordOffs[data[0][0]]=[];
		if (dlength<0) {
			//For deletions relative to reference
			while (coordOffs[data[0][0]].length<data[0][1][1]) coordOffs[data[0][0]].push(0);
			for (var i=data[0][1][0]; i<data[0][1][1]; i++) coordOffs[data[0][0]][i]=-1;
			while (coordOffs[data[0][0]].length<data[0][1][1]) coordOffs[data[0][0]].push(-1);
		} else if (dlength>0) {
			//For insertions relative to reference
			while (coordOffs[data[0][0]].length<data[0][1][0]) coordOffs[data[0][0]].push(0);
			coordOffs[data[0][0]][data[0][1][1]]=dlength;
		}
	}
});

//On end of input
rl.on("close", function() {
	//Determine cumulative coordinate remap
	for (contig in coordOffs) {
		coordMaps[contig]=[];
		var index=0;
		for (var i=0; i<coordOffs[contig].length; i++) {
			coordMaps[contig][i]=index;
			index+=1+coordOffs[contig][i];
		}
	}
	//Remap GFF coordinates and print
	for (var i=0; i<gff.length; i++) {
		if (gff[i].length>2 && gff[i].substring(0, 2)!="##") {
			gff[i]=gff[i].split("\t");
			if (coordMaps[gff[i][0]]) {
				//If this contig is remapped (otherwise do not modify coordinates)
				for (var j=0; j<2; j++) {
					var coord=parseInt(gff[i][3+j]);
					//Extend coordinate mapping if necessary
					while (coord>=coordMaps[gff[i][0]].length) coordMaps[gff[i][0]].push(coordMaps[gff[i][0]][coordMaps[gff[i][0]].length-1]+1);
					gff[i][3+j]=coordMaps[gff[i][0]][coord];
				}
			}
			console.log(gff[i].join("\t"));
		}
	}
});