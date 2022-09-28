var fs=require("fs");
var path=require("path");

var type="nucl";
if (process.argv[3]=="prot") type="prot";
console.error("Outputting "+type+" sequences");

//https://doi.org/10.1016/0014-4894(91)90012-L
var codon={
	"A": [
		{"prob": 0.27, "codon": "GCT"},
		{"prob": 0.25, "codon": "GCC"},
		{"prob": 0.25, "codon": "GCA"},
		{"prob": 0.25, "codon": "GCG"}
	],
	"R": [
		{"prob": 0.27, "codon": "CGT"},
		{"prob": 0.31, "codon": "CGC"},
		{"prob": 0.11, "codon": "CGA"},
		{"prob": 0.14, "codon": "CGG"},
		{"prob": 0.06, "codon": "AGA"},
		{"prob": 0.11, "codon": "AGG"}
	],
	"N": [
		{"prob": 0.39, "codon": "AAT"},
		{"prob": 0.61, "codon": "AAC"}
	],
	"D": [
		{"prob": 0.48, "codon": "GAT"},
		{"prob": 0.52, "codon": "GAC"}
	],
	"C": [
		{"prob": 0.36, "codon": "TGT"},
		{"prob": 0.64, "codon": "TGC"}
	],
	"Q": [
		{"prob": 0.37, "codon": "CAA"},
		{"prob": 0.63, "codon": "CAG"}
	],
	"E": [
		{"prob": 0.34, "codon": "GAA"},
		{"prob": 0.66, "codon": "GAG"}
	],
	"G": [
		{"prob": 0.42, "codon": "GGT"},
		{"prob": 0.28, "codon": "GGC"},
		{"prob": 0.17, "codon": "GGA"},
		{"prob": 0.14, "codon": "GGG"}
	],
	"H": [
		{"prob": 0.36, "codon": "CAT"},
		{"prob": 0.64, "codon": "CAC"}
	],
	"I": [
		{"prob": 0.44, "codon": "ATT"},
		{"prob": 0.41, "codon": "ATC"},
		{"prob": 0.15, "codon": "ATA"}
	],
	"L": [
		{"prob": 0.06, "codon": "TTA"},
		{"prob": 0.15, "codon": "TTG"},
		{"prob": 0.26, "codon": "CTT"},
		{"prob": 0.22, "codon": "CTC"},
		{"prob": 0.08, "codon": "CTA"},
		{"prob": 0.23, "codon": "CTG"}
	],
	"K": [
		{"prob": 0.30, "codon": "AAA"},
		{"prob": 0.70, "codon": "AAG"}
	],
	"F": [
		{"prob": 0.45, "codon": "TTT"},
		{"prob": 0.55, "codon": "TTC"}
	],
	"P": [
		{"prob": 0.22, "codon": "CCT"},
		{"prob": 0.29, "codon": "CCC"},
		{"prob": 0.28, "codon": "CCA"},
		{"prob": 0.22, "codon": "CCG"}
	],
	"S": [
		{"prob": 0.17, "codon": "TCT"},
		{"prob": 0.21, "codon": "TCC"},
		{"prob": 0.16, "codon": "TCA"},
		{"prob": 0.14, "codon": "TCG"},
		{"prob": 0.16, "codon": "AGT"},
		{"prob": 0.16, "codon": "AGC"}
	],
	"T": [
		{"prob": 0.21, "codon": "ACT"},
		{"prob": 0.25, "codon": "ACC"},
		{"prob": 0.25, "codon": "ACA"},
		{"prob": 0.29, "codon": "ACG"}
	],
	"Y": [
		{"prob": 0.35, "codon": "TAT"},
		{"prob": 0.65, "codon": "TAC"}
	],
	"V": [
		{"prob": 0.29, "codon": "GTT"},
		{"prob": 0.16, "codon": "GTC"},
		{"prob": 0.14, "codon": "GTA"},
		{"prob": 0.41, "codon": "GTG"}
	],
	"M": [
		{"prob": 1.00, "codon": "ATG"}
	],
	"W": [
		{"prob": 1.00, "codon": "TGG"}
	],
	"*": [
		{"prob": 0.48, "codon": "TAA"},
		{"prob": 0.26, "codon": "TAG"},
		{"prob": 0.26, "codon": "TGA"}
	]
}

var invcodon={};
for (var key in codon) {
	for (var i=0; i<codon[key].length; i++) {
		invcodon[codon[key][i].codon]=key;
	}
}

//Load GFF file
filepath=path.resolve(process.cwd(), process.argv[2]);
var data=fs.readFileSync(filepath).toString().replace(/\r/g, "");
data=data.split("\n");
//Loop through lines to find fasta data
contigs=[];
for (var i=0; i<data.length; i++) {
	if (data[i][0]==">") {
		var name=data[i].substring(1);
		i++;
		contigs[name]=data[i];
	}
}
//Loop through lines to find CDSs to extract
for (var i=0; i<data.length; i++) {
	if (data[i][0]!="#" && data[i][0]!=">") {
		data[i]=data[i].split("\t");
		if (data[i][2]=="CDS") {
			data[i][8]=data[i][8].split(";");
			var metaObj={};
			for (var j in data[i][8]) {
				data[i][8][j]=data[i][8][j].split("=");
				metaObj[data[i][8][j][0]]=data[i][8][j][1];
			}
			data[i][8]=metaObj;
			console.log([">"+data[i][8].Parent, "|", data[i][0], data[i][3], data[i][4], data[i][6]].join(" "));
			if (data[i][6]=="+") {
				var seq=contigs[data[i][0]].substring(parseInt(data[i][3])-1, parseInt(data[i][4]));
			} else if (data[i][6]=="-") {
				var seq=contigs[data[i][0]].substring(parseInt(data[i][3])-1, parseInt(data[i][4]));
				seq=reverseComplement(seq);
			}
			if (type=="prot") seq=translate(seq);
			console.log(seq);
		}
	} else if (data[i][0]==">") {
		//To skip fasta data
		i=data.length;
	}
}

function translate(seq) {
	var oseq="";
	for (var i=0; i<seq.length; i+=3) {
		oseq+=invcodon[seq.substring(i, i+3)];
	}
	return oseq;
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