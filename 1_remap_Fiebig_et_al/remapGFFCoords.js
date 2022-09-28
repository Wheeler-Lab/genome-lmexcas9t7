var fs=require("fs");
var path=require("path");

remap=require(path.resolve(process.cwd(), process.argv[3]));

var gff=(fs.readFileSync(path.resolve(process.cwd(), process.argv[2])).toString()).replace(/\r/g, "");
var gff=gff.split("\n");
for (var i=0; i<gff.length; i++) {
	if (gff[i].length==0) {
		console.log(gff[i]);
	} else {
		if (gff[i].substring(0, "##".length)=="##") {
			console.log(gff[i]);
		} else {
			gff[i]=gff[i].split("\t");
			gff[i][3]=parseInt(gff[i][3]);
			gff[i][4]=parseInt(gff[i][4]);
			//Check for remapping
			for (var contig in remap) {
				if (remap[contig].offs!=0) {
					//Only check if this contig's offset is not zero
					if (gff[i][0]==remap[contig].contig && gff[i][3]>=remap[contig].start+remap[contig].offs && gff[i][3]<=remap[contig].end+remap[contig].offs && gff[i][4]>=remap[contig].start+remap[contig].offs && gff[i][4]<=remap[contig].end+remap[contig].offs) {
						//If the annotation sits fully within the remapped range
						gff[i][0]=contig;
						gff[i][3]-=remap[contig].offs;
						gff[i][4]-=remap[contig].offs;
					}
				}
			}
			console.log(gff[i].join("\t"));
		}
	}
}