var fs=require("fs");
var path=require("path");

var gff=(fs.readFileSync(path.resolve(process.cwd(), process.argv[2])).toString()).replace(/\r/g, "");
gff=gff.split("\n");
var ids=[];
//For each line
for (var i=1; i<gff.length; i++) {
	//Split to columns
	gff[i]=gff[i].split(",");
	if (gff[i].length>=8) {
		gff[i][2]="gene";
		//Do not add duplicate IDs
		if (ids.indexOf(gff[i][8])==-1) {
			ids.push(gff[i][8]);
			//Set up additional parameters (col 8)
			gff[i][8]={
				ID: gff[i][8]
			};
			//Print result
			//gene, mRNA (rna_), CDS (cds_), exon (exon_)
			var params=[];
			for (var param in gff[i][8]) params.push(""+param+"="+gff[i][8][param]);
			var cur=gff[i].slice(0, 8);
			cur[cur.length]=params.join(";");
			console.log(cur.join("\t"));
			var sub={
				"mRNA": {
					prefix: "rna_",
					suffix: "-1",
					parentPrefix: "",
					parentSuffix: ""
				},
				"CDS": {
					prefix: "cds_",
					suffix: "-1",
					parentPrefix: "rna_",
					parentSuffix: "-1"
				},
				"exon": {
					prefix: "exon_",
					suffix: "-1",
					parentPrefix: "rna_",
					parentSuffix: "-1"
				}
			};
			for (s in sub) {
				var cur=gff[i].slice(0, 8);
				cur[2]=s;
				if (s=="CDS") cur[7]=0;
				var curparams={};
				for (param in gff[i][8]) curparams[param]=gff[i][8][param];
				curparams.ID=sub[s].prefix+gff[i][8].ID+sub[s].suffix;
				curparams.Parent=sub[s].parentPrefix+gff[i][8].ID+sub[s].parentSuffix;
				var curparamsa=[];
				for (var param in curparams) curparamsa.push(""+param+"="+curparams[param]);
				cur[cur.length]=curparamsa.join(";");
				console.log(cur.join("\t"));
			}
		}
	}
}