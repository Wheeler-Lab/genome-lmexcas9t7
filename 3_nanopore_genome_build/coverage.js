var curid="";
var curcov=0;
var curlen=0;
var readline=require("readline");
var rl=readline.createInterface({
	input: process.stdin,
	output: process.stdout,
	terminal: false
});

rl.on("line", function(line){
	if (curid!=line.split("\t")[0]) {
		if (curid!="") {
			console.log(curid, curcov/curlen, curlen);
		}
		curid=line.split("\t")[0];
		curcov=0;
		curlen=0;
	}
	curcov+=parseInt(line.split("\t")[2]);
	curlen=parseInt(line.split("\t")[1]);
});

rl.on("close", function() {
	if (curid!="") {
		console.log(curid, curcov/curlen, curlen);
	}
});