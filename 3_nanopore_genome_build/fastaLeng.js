var readline=require("readline");

var id="";
var len="";

var rl=readline.createInterface({
	input: process.stdin,
	output: process.stdout,
	terminal: false
});

rl.on("line", function(line){
	if (line[0]==">") {
		if (id!="") {
			console.log(id, len);
		}
		id=line.substring(1, line.length);
		len=0;
	} else {
		len+=line.length;
	}
});

rl.on("close", function() {
	console.log(id, len);
});