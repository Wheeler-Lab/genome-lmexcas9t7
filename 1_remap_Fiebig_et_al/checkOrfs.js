var fs=require("fs");
var path=require("path");
var codonl=3;

//Load fasta file
filepath=path.resolve(process.cwd(), process.argv[3]);
var fasta=("\n"+(fs.readFileSync(filepath).toString()).replace(/\r/g, "")).split(">");
contigs=[];
for (var i=1; i<fasta.length; i++) {
	var name=fasta[i].substring(0, fasta[i].indexOf("\n"));
	if (name.indexOf(" ")!=-1) name=name.substring(0, name.indexOf(" "));
	var seq=fasta[i].substring(fasta[i].indexOf("\n")).replace(/\n/gi, "");
	contigs[name]=seq;
}

//ORF-checking functions
function checkOrf(seq) {
	var start=["ATG"];
	var stop=["TAG", "TAA", "TGA"];
	var res={
		start: true,
		stop: true,
		orf: true,
		ns: false
	};
	if (start.indexOf(seq.substring(0, codonl))==-1) res.start=false;
	if (stop.indexOf(seq.substring(seq.length-codonl))==-1) res.stop=false;
	if (seq.length%codonl!=0) res.orf=false;
	for (var i=0; i<seq.length-codonl; i+=codonl) if (stop.indexOf(seq.substring(i, i+codonl))!=-1) res.orf=false;
	if (seq.indexOf("N")!=-1) res.ns=true;
	return res;
}

//Fetch sequence from contigs
function fetchSequence(start, stop, dir, contig, contigs) {
	if (contigs[contig]) {
		var seq=contigs[contig].substring(Math.max(parseInt(start)-1, 0), Math.min(parseInt(stop), contigs[contig].length)).toUpperCase();
		if (dir=="-") seq=reverseComplement(seq);
		return seq;
	} else {
		console.error(contig+" not found");
		return "";
	}
}

//Reverse complement IUPAC sequence
function reverseComplement(seq) {
	var f=["A", "C", "G", "T", "U", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V", "N", ".", "-"];
	var r=["T", "G", "C", "A", "A", "Y", "R", "S", "W", "M", "K", "V", "H", "D", "B", "N", ".", "-"];
	var oseq="";
	for (var i=0; i<seq.length; i++) {
		var cur=seq.substring(i, i+1).toUpperCase();
		if (f.indexOf(cur)!=-1) {
			oseq=r[f.indexOf(cur)]+oseq;
		}
	}
	return oseq;
}

//Helper function to find children of specified type
function findChildrenIndex(arr, type) {
	var index=[];
	for (var i in arr) if (arr[i].type==type) index.push(i);
	return index;
}

//Load JSON file
json=require(path.resolve(process.cwd(), process.argv[2]));

for (entry in json) {
	//For each json entry
	if (json[entry].type=="gene") {
		//If it is a gene
		var mrnaIndex=findChildrenIndex(json[entry].children, "mRNA");
		for (var i in mrnaIndex) {
			//and it has an mRNA child
			var cdsIndex=findChildrenIndex(json[entry].children[mrnaIndex[i]].children, "CDS");
			//Determine the coordinate extent of the CDSs
			var start=false;
			var end=false;
			for (var j in cdsIndex) {
				if (!start) start=json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].start;
				if (!end) end=json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].end;
				start=Math.min(start, json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].start);
				end=Math.min(end, json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].end);
				//Check each CDS individually while determining coordinate extent of CDSs
				var seq=fetchSequence(json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].start, json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].end, json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].dir, json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].contig, contigs);
				var res=checkOrf(seq);
				json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].params.startCodon=res.start;
				json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].params.stopCodon=res.stop;
				json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].params.clearORF=res.orf;
				json[entry].children[mrnaIndex[i]].children[cdsIndex[j]].params.ambiguousSequence=res.ns;
			}
			//If a coordinate extent can be determined
			if (start && end) {
				//Fetch the sequence
				//NB: Currently collapses all CDSs for an mRNA to a single sequence
				//Ignores introns, ignores posibility of multiple CDSs per mRNA
				var seq=fetchSequence(start, end, json[entry].children[mrnaIndex[i]].dir, json[entry].children[mrnaIndex[i]].contig, contigs);
				//Check the sequence and sets the gene parameters
				var res=checkOrf(seq);
				json[entry].params.startCodon=res.start;
				json[entry].params.stopCodon=res.stop;
				json[entry].params.clearORF=res.orf;
				json[entry].params.ambiguousSequence=res.ns;
			}
		}
	}
}

console.log(JSON.stringify(json, "", "\t"));