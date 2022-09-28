var fs=require("fs");
var path=require("path");

//Load GFF file
filepath=path.resolve(process.cwd(), process.argv[2]);
var gff=fs.readFileSync(filepath).toString().replace(/\r/g, "");
gff=gff.split("\n");

//Load fasta file
filepath=path.resolve(process.cwd(), process.argv[3]);
var fasta="\n"+(fs.readFileSync(filepath).toString()).replace(/\r/g, "");
fasta=fasta.split(">");
contigs=[];
for (var i=1; i<fasta.length; i++) {
	var name=fasta[i].substring(0, fasta[i].indexOf("\n"));
	if (name.indexOf(" ")!=-1) name=name.substring(0, name.indexOf(" "));
	var seq=fasta[i].substring(fasta[i].indexOf("\n")).replace(/\n/gi, "");
	contigs[name]=seq;
}

//Print output
console.log(["##gff-version", 3].join("\t"));
for (var contig in contigs) console.log(["##sequence-region", contig, 1, contigs[contig].length].join("\t"));
for (var line in gff) console.log(gff[line]);
console.log("##FASTA")
for (var contig in contigs) {
	console.log(">"+contig);
	console.log(contigs[contig]);
}