var readline=require("readline");

var histo=[];
var type="";
var running=0;
var index=0;

var rl=readline.createInterface({
	input: process.stdin,
	output: process.stdout,
	terminal: false
});

rl.on("line", function(line){
	if (type=="" && line.length>1) {
		type=line[0];
		if (type=="@") {
			index=0;
		}
		if (type!="@" && type!=">") {
			console.error("Warning: Unclear if input is FASTA or FASTQ");
			console.error("         Assuming FASTA and continuing");
			type=">";
		}
	}
	if (type==">") {
		if (line[0]==type) {
			if (running!=0) {
				for (var i=histo.length; i<=running; i++) {
					histo[i]=0;
				}
				histo[running]+=1;
			}
			running=0;
		} else {
			running+=line.length;
		}
	} else if (type=="@") {
		if (index%4==1) {
			for (var i=histo.length; i<=line.length; i++) {
				histo[i]=0;
			}
			histo[line.length]+=1;
			pnext=false;
		}
		index++;
	}
});

rl.on("close", function() {
	for (var i=0; i<histo.length; i++) {
		console.log(i, histo[i], i*histo[i]);
	}
});