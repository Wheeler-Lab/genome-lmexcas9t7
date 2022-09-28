var fs=require("fs");
var path=require("path");

function loadFasta(file) {
	filepath=path.resolve(process.cwd(), file);
	var fasta=("\n"+(fs.readFileSync(filepath).toString()).replace(/\r/g, "")).split(">");
	contigs=[];
	for (var i=1; i<fasta.length; i++) {
		var name=fasta[i].substring(0, fasta[i].indexOf("\n"));
		if (name.indexOf(" ")!=-1) name=name.substring(0, name.indexOf(" "));
		var seq=fasta[i].substring(fasta[i].indexOf("\n")).replace(/\n/gi, "");
		contigs[name]=seq;
	}
	return contigs;
}

//Original JSON
json=require(path.resolve(process.cwd(), process.argv[2]));
//Additional JSON
addjson=require(path.resolve(process.cwd(), process.argv[3]));
//Reference FASTA
fasta=loadFasta(path.resolve(process.cwd(), process.argv[4]));

function checkOrf(seq) {
	var codonl=3;
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

function fetchCdsSequence(mrna, contigs) {
	mrna.children.sort();
	mrna.children.sort(function(a, b) {
			return a.start-b.start;
	});
	seq="";
	for (child in mrna.children) if (mrna.children[child].type=="CDS") seq+=contigs[mrna.children[child].contig].substring(mrna.children[child].start-1, mrna.children[child].end).toUpperCase();
	if (mrna.children[child].dir=="-") seq=reverseComplement(seq);
	return seq;
}

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

function getTypeIndex(arr, type) {
	for (var i=0; i<arr.length; i++) if (arr[i].type==type) return i;
	return -1;
}

for (var addentry in addjson) {
	var entryindex=-1;
	for (var entry in json) if (json[entry].id==addjson[addentry].id) entryindex=entry;
	if (entryindex!=-1) {
		//Fetch old and new CDSs and check sequences
		if (!addjson[addentry].children) console.error(addjson[addentry]);
		//Only considers the first mrna entry
		var mrna=getTypeIndex(addjson[addentry].children, "mRNA");
		var seq=fetchCdsSequence(addjson[addentry].children[mrna], fasta);
		var res=checkOrf(seq);
		if (res.start && res.stop && res.orf & json[entryindex].dir==addjson[addentry].dir) {
			//Replace mRNA entry if an orf and on the same strand
			console.error(json[entryindex].id+" CDS updated");
			json[entryindex].params.cdsHistory="Updated_FiebigEtAl2015";
			var rmrna=getTypeIndex(json[entryindex].children, "mRNA");
			json[entryindex].children[rmrna]=addjson[addentry].children[mrna];
		} else {
			console.error(json[entryindex].id+" CDS not updated");
		}
	}
}

console.log(JSON.stringify(json, "", "\t"));