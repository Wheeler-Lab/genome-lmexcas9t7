var readline=require("readline");
var snp=0;
var ins=0;
var del=0;
var oth=0;

var rl=readline.createInterface({
	input: process.stdin,
	output: process.stdout,
	terminal: false
});

rl.on("line", function(line){
	var data=line.split(" ");
	console.error(line);
	if (data[2]==".") {
		ins++;
	} else if (data[3]==".") {
		del++;
	} else if (data[2].length==1 && data[3].length==1) {
		snp++;
	} else if (data[2]!=data[3]) {
		oth++;
	}
});

rl.on("close", function() {
	console.log(snp, ins, del, oth);
});