var fs=require("fs");
var path=require("path");

function loadGFF(file) {
	filepath=path.resolve(process.cwd(), file);
	var gff=(fs.readFileSync(filepath).toString()).replace(/\r/g, "");
	gff=gff.split("\n");
	for (var i=0; i<gff.length; i++) {
		if (gff[i].length==0) {
			gff.splice(i, 1);
			i--;
		} else if (gff[i].substring(0, "#".length)=="#") {
			gff.splice(i, 1);
			i--;
		} else if (gff[i][0]==">") {
			gff.splice(i, gff.length-i);
		} else {
			gff[i]=gff[i].split("\t");
			gff[i][3]=parseInt(gff[i][3]);
			gff[i][4]=parseInt(gff[i][4]);
			gff[i][8]=gff[i][8].split(";");
			var params={};
			for (var j=0; j<gff[i][8].length; j++) {
				gff[i][8][j]=gff[i][8][j].split("=");
				params[gff[i][8][j][0]]=gff[i][8][j][1];
			}
			gff[i][8]=params;
		}
	}
	return gff;
}

var gff=loadGFF(path.resolve(process.cwd(), process.argv[2]));
var json=[];

//Function for recursively finding objects by ID and adding as child
function recursiveAddToJson(json, parentId, res) {
	for (entry in json) {
		if (json[entry].id==parentId) {
			if (!json[entry].children) json[entry].children=[];
			json[entry].children.push(res);
			return true;
		}
		if (json[entry].children) {
			var success=recursiveAddToJson(json[entry].children, parentId, res);
			if (success) return true;
		}
	}
	return false;
}

//Parse GFF to json
for (var i=0; i<gff.length; i++) {
	//Only add the object to the json if it has an ID
	if (gff[i][8].ID) {
		//Only continue if this ID has not already been used
		//Make an object recording the gff entry
		var res={
			id: gff[i][8].ID,
			type: gff[i][2],
			contig: gff[i][0],
			start: gff[i][3],
			end: gff[i][4],
			dir: gff[i][6],
			params: gff[i][8]
		};
		if (gff[i][7]!=".") res.frame=parseInt(gff[i][7]);
		if (gff[i][8].Parent) {
			//If it has a parent then find and add as a child to the parent
			success=recursiveAddToJson(json, gff[i][8].Parent, res);
			if (!success) {
				//Assume the parent is a missing root level element
				console.error("Parent ID "+gff[i][8].Parent+" for "+gff[i][8].ID+" not found");
				json.push({
					id: gff[i][8].Parent,
					children: [res]
				});
			}
		} else {
			//Otherwise add as a new json entry
			json.push(res);
		}
	}
}
console.log(JSON.stringify(json, "", "\t"));