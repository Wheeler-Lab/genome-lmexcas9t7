var fs=require("fs");
var path=require("path");

var gff=(fs.readFileSync(path.resolve(process.cwd(), process.argv[2])).toString()).replace(/\r/g, "");
gff=gff.split("\n");
//For each line
for (var i=0; i<gff.length; i++) {
	//Split to columns
	gff[i]=gff[i].split(",");
	if (gff[i].length>=8) {
		//Split additional parameters (col 8)
		gff[i][8]=gff[i][8].split(";");
		var params={};
		for (var j=0; j<gff[i][8].length; j++) {
			if (gff[i][8][j].length==0) {
				gff[i][8].splice(j, 1);
				j--;
			} else {
				gff[i][8][j]=gff[i][8][j].split("=");
				params[gff[i][8][j][0]]=gff[i][8][j][1];
			}
		}
		gff[i][8]=params;
		//Derive Parent from ID parameter
		var pid=gff[i][8].ID.split("_");
		pid.splice(pid.length-2, 2);
		pid=pid.join("_");
		gff[i][8].Parent=pid;
		//Print result
		var params=[];
		for (var param in gff[i][8]) params.push(""+param+"="+gff[i][8][param]);
		gff[i][8]=params.join(";");
		console.log(gff[i].join("\t"));
	}
}