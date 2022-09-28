var fs=require("fs");
var path=require("path");

function loadGFF(file, out) {
	var gff=(fs.readFileSync(path.resolve(process.cwd(), file)).toString()).replace(/\r/g, "");
	var gff=gff.split("\n");
	var seqs={};
	for (var i=0; i<gff.length; i++) {
		if (gff[i].length==0) {
			gff.splice(i, 1);
			i--;
		} else if (gff[i].substring(0, "#".length)=="#") {
			if (gff[i].indexOf("\t")!=-1) {
				gff[i]=gff[i].split("\t");
			} else {
				gff[i]=gff[i].split(" ");
			}
			if (gff[i][0]=="##sequence-region") {
				seqs[gff[i][1]]={
					start: parseInt(gff[i][2]),
					end: parseInt(gff[i][3])
				};
			}
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
			//if (!gff[i][8].ID) {
			//	gff.splice(i, 1);
			//	i--;
			//}
		}
	}
	if (out=="seqs") {
		return seqs;
	}
	return gff;
}

var gffold=loadGFF(process.argv[2], "gff");
var seqsold=loadGFF(process.argv[2], "seqs");
var gffref=loadGFF(process.argv[3], "gff");
var seqsref=loadGFF(process.argv[3], "seqs");

//Determine coordinate remapping for sequences
var remap={};
for (var contig in seqsref) {
	//For each contig in the reference annotation
	var refgenes=[];
	for (var generef in gffref) if (gffref[generef][0]==contig) refgenes.push(generef);
	//Find a gene matching the ID of one of the reference genes in the old annotation
	var oldindex=-1;
	var refindex=-1;
	var refgeneindex=0;
	for (refgeneindex=0; refgeneindex<refgenes.length; refgeneindex++) {
		for (var geneold in gffold) if (gffold[geneold][8].ID==gffref[refgenes[refgeneindex]][8].ID) {
			oldindex=geneold;
			refindex=refgenes[refgeneindex];
		}
	}
	if (refindex!=-1 && oldindex!=-1) {
		var offs=gffold[oldindex][3]-gffref[refindex][3];
		console.error(contig+" remapped from "+gffold[oldindex][0]+" with "+offs+" base offset (based on "+gffref[refindex][8].ID+")");
		remap[contig]={
			contig: gffold[oldindex][0],
			offs: offs,
			start: seqsref[contig].start+offs,
			end: seqsref[contig].end+offs
		};
	} else {
		console.error(contig+" not remapped (no matching geneIDs on this contig)");
	}
}

console.log(JSON.stringify(remap, "", "\t"));