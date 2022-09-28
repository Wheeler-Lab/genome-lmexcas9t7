var fs=require("fs");
var path=require("path");

var offsLimitStart=80;
var offsLimitStop=200;
var codonl=3;

var verbose=true;
var indent="\t";

//Load GFF file
filepath=path.resolve(process.cwd(), process.argv[2]);
var gff=fs.readFileSync(filepath).toString().replace(/\r/g, "").split("\n");

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

//Loop through GFF lines
var i=0;
//Until a FASTA data block is hit
while (gff[i][0]!=">" && i<gff.length-1) {
	//If not a comment
	var comment="##";
	if (gff[i].substring(0, comment.length)!=comment) {
		gff[i]=gff[i].split("\t");
		//Match criteria for which features to correct
		if (gff[i].length>=8) {
			var id=gff[i][2];
			gff[i][3]=parseInt(gff[i][3]);
			gff[i][4]=parseInt(gff[i][4]);
			var seq=fetchSequence(gff[i][3], gff[i][4], gff[i][6], gff[i][0]);
			var res=checkOrf(seq);
			if (res.ns && verbose) console.error(indent+"Warning: "+id+" has Ns");
			if (res.start && res.stop && res.orf) {
				//Do not modify genes with start, stop and orf
				console.log(gff[i].join("\t"));
			} else {
				if (seq.length%codonl!=0) {
					//Check for non-integer number of codons
					if (verbose) console.error(indent+"Warning: "+id+" trimmed to integer number of codons");
					var length=gff[i][4]-gff[i][3];
					length=Math.floor(length/codonl)*codonl;
					if (res.start) {
						//Trim end if start codon present
						if (gff[i][6]=="+") {
							gff[i][4]=gff[i][3]+length;
						} else if (gff[i][6]=="-") {
							gff[i][3]=gff[i][4]-length;
						}
					} else {
						//Otherwise trim start
						if (gff[i][6]=="+") {
							gff[i][3]=gff[i][4]-length;
						} else if (gff[i][6]=="-") {
							gff[i][4]=gff[i][3]+length;
						}
					}
				}
				//Attempt to fix genes with an integer number of codons
				if (!res.orf) {
					//If orf has stop codons in it
					//Attempt to trim to premature stop codons
					offs=-1;
					while (!res.orf && offs*codonl<seq.length && offs<offsLimitStop) {
						offs++
						if (gff[i][6]=="+") {
							seq=fetchSequence(gff[i][3], gff[i][4]-offs*codonl, gff[i][6], gff[i][0]);
						} else if (gff[i][6]="-") {
							seq=fetchSequence(gff[i][3]+offs*codonl, gff[i][4], gff[i][6], gff[i][0]);
						}
						res=checkOrf(seq);
					}
					if (res.orf) {
						if (gff[i][6]=="+") {
							console.error(id+" stop ("+gff[i][0]+" "+gff[i][4]+" "+gff[i][6]+") trimmed by "+offs+" codons");
							gff[i][4]=gff[i][4]-offs*codonl;
						} else if (gff[i][6]=="-") {
							console.error(id+" stop ("+gff[i][0]+" "+gff[i][3]+" "+gff[i][6]+") trimmed by "+offs+" codons");
							gff[i][3]=gff[i][3]+offs*codonl;
						}
					} else {
						if (verbose) console.error(indent+"Error: "+id+" premature stop codons could not be corrected");
					}
				} else if (!res.stop) {
					//Otherwise, if a stop codon is missing
					//Attempt to extend to find a stop codon
					if (!res.stop) {
						offs=-1;
						while (!res.stop && offs<offsLimitStop) {
							offs++;
							if (gff[i][6]=="+") {
								seq=fetchSequence(gff[i][3], gff[i][4]+offs*codonl, gff[i][6], gff[i][0]);
							} else if (gff[i][6]=="-") {
								seq=fetchSequence(gff[i][3]-offs*codonl, gff[i][4], gff[i][6], gff[i][0]);
							}
							res=checkOrf(seq);
						}
						if (res.stop && res.orf) {
							if (gff[i][6]=="+") {
								console.error(id+" stop ("+gff[i][0]+" "+gff[i][4]+" "+gff[i][6]+") extended by "+offs+" codons");
								gff[i][4]=gff[i][4]+offs*codonl;
							} else if (gff[i][6]=="-") {
								console.error(id+" stop ("+gff[i][0]+" "+gff[i][3]+" "+gff[i][6]+") extended by "+offs+" codons");
								gff[i][3]=gff[i][3]-offs*codonl;
							}
						} else {
							if (verbose) console.error(indent+"Error: "+id+" stop codon could not be corrected");
						}
					}
				}
				if (!res.start) {
					//Attempt to fix start
					//Scan for upstream starts
					uoffs=-1;
					ures=checkOrf(seq);
					while (!ures.start && uoffs<offsLimitStart) {
						uoffs++;
						if (gff[i][6]=="+") {
							useq=fetchSequence(gff[i][3]-uoffs*codonl, gff[i][4], gff[i][6], gff[i][0]);
						} else if (gff[i][6]=="-") {
							useq=fetchSequence(gff[i][3], gff[i][4]+uoffs*codonl, gff[i][6], gff[i][0]);
						}
						ures=checkOrf(useq);
					}
					//Scan for downstream starts
					doffs=-1;
					dres=checkOrf(seq);
					while (!dres.start && doffs*codonl<seq.length && doffs<offsLimitStart) {
						doffs++;
						if (gff[i][6]=="+") {
							dseq=fetchSequence(gff[i][3]+doffs*codonl, gff[i][4], gff[i][6], gff[i][0]);
						} else if (gff[i][6]=="-") {
							dseq=fetchSequence(gff[i][3], gff[i][4]-doffs*codonl, gff[i][6], gff[i][0]);
						}
						dres=checkOrf(dseq);
					}
					//Decide whether to use up or downstream starts
					var offs=uoffs;
					var res=ures;
					var mod="extended";
					if (uoffs>doffs) {
						offs=-doffs;
						res=dres;
						mod="trimmed";
					}			
					if (res.start && res.orf) {
						if (gff[i][6]=="+") {
							console.error(id+" start ("+gff[i][0]+" "+gff[i][3]+" "+gff[i][6]+") "+mod+" by "+Math.abs(offs)+" codons");
							gff[i][3]=gff[i][3]-offs*codonl;
						} else if (gff[i][6]=="-") {
							console.error(id+" start ("+gff[i][0]+" "+gff[i][4]+" "+gff[i][6]+") "+mod+" by "+Math.abs(offs)+" codons");
							gff[i][4]=gff[i][4]+offs*codonl;
						}
					} else {
						if (verbose) console.error(indent+"Error: "+id+" start codon could not be corrected");
					}
				}
				//Print the corrected coordinates
				console.log(gff[i].join("\t"));
			}
		} else {
			//Print uncorrected features
			console.log(gff[i].join("/t"));
		}
	} else {
		//Print comments
		console.log(gff[i]);
	}
	i++;
}
//Roll back index one then print fasta data block
i--;
while (i<gff.length) {
	console.log(gff[i]);
	i++;
}

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

function fetchSequence(start, stop, dir, contig) {
	var seq=contigs[contig].substring(Math.max(parseInt(start)-1, 0), Math.min(parseInt(stop), contigs[contig].length)).toUpperCase();
	if (dir=="-") seq=reverseComplement(seq);
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