var fs=require("fs");
var path=require("path");

json=require(path.resolve(process.cwd(), process.argv[2]));

function recursivePrintGff(json) {
	for (var i=0; i<json.length; i++) {
		var params=[];
		for (var entry in json[i].params) params.push(""+entry+"="+json[i].params[entry]);
		if (!("frame" in json[i]))
			json[i].frame=".";
		console.log([json[i].contig, "RW1", json[i].type, json[i].start, json[i].end, ".", json[i].dir, json[i].frame, params.join(";")].join("\t"));
		if (json[i].children) recursivePrintGff(json[i].children);
	}
}

recursivePrintGff(json);